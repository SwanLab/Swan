classdef DifferentSubdomainsEIFEM < handle

    properties (Access = public)

    end
    properties (Access = private)
        referenceMesh
        meshDomain
        boundaryConditions
        bcApplier
        LHS
        RHS
        fileNameEIFEM
        tolSameNode
        solverType
    end


    properties (Access = private)
        nSubdomains
    end

    methods (Access = public)

        function obj = DifferentSubdomainsEIFEM()
            close all
            obj.init()
            %radius_to_analyse = 0.05:0.005:0.7;
            %Kc = cell(length(radius_to_analyse),1);
            radiusMesh = 0.5;

            obj.createReferenceMesh(false,radiusMesh);
            bS  = obj.referenceMesh.createBoundaryMesh();
            [mD,mSb,iC,lG,iCR,discMesh] = obj.createMeshDomain();
            obj.meshDomain = mD;
            [bC,dir] = obj.createBoundaryConditions();
            obj.boundaryConditions = bC;
            obj.createBCapplier()

            [LHSr, Mr, RHSr] = obj.createElasticProblem();

                        
            LHSfun = @(x) LHSr*x;
            [Meifem, Kcoarse, Mcoarse, EIFEM, ss]       = obj.createEIFEMPreconditioner(dir,iC,lG,bS,iCR,discMesh,radiusMesh);
            
            Milu         = obj.createILUpreconditioner(LHSr);
            Mmult        = @(r) Preconditioner.multiplePrec(r,LHSfun,Milu,Meifem,Milu);
            Mid          = @(r) r;

            tol = 1e-8;
            x0 = zeros(size(RHSr));
            xSol = LHSr\RHSr;

            [uPCG,residualPCG,errPCG,errAnormPCG] = PCG.solve(LHSfun,RHSr,x0,Mmult,tol,xSol);     
            [uCG,residualCG,errCG,errAnormCG]    = PCG.solve(LHSfun,RHSr,x0,Mid,tol,xSol);

            uCGFull = obj.bcApplier.reducedToFullVectorDirichlet(uCG); %reapply all dirichlett DOFs
            uF = saveDeformed(obj.meshDomain,uCGFull);
            plot(uF)
        end

    end

    methods (Access = private)

        function init(obj)
            obj.nSubdomains  = [15 5]; %nx ny
            obj.fileNameEIFEM = 'DEF_Q4porL_1.mat';
            obj.tolSameNode = 1e-10;
            obj.solverType = 'REDUCED';
        end        

        function createReferenceMesh(obj, loadfile, radius)
            if loadfile
                filename = obj.fileNameEIFEM;
                load(filename); 
                s.coord     = EIFEoper.MESH.COOR;
                s.connec    = EIFEoper.MESH.CN;
                s.interType = 'QUADRATIC';
                
            else
                holeMesh    = obj.createMesh(radius);      
                s.coord     = holeMesh.coord;
                s.connec    = holeMesh.connec;
                s.interType = 'LINEAR';
                s           = obj.updateCoordsMesh(s); 
            end
        
            obj.referenceMesh = Mesh.create(s);
        end
        
        function holeMesh = createMesh(obj,radius)
            fullmesh = UnitTriangleMesh(12,12);
            ls = obj.computeCircleLevelSet(fullmesh,radius);
            sUm.backgroundMesh = fullmesh;
            sUm.boundaryMesh   = fullmesh.createBoundaryMesh;
            uMesh              = UnfittedMesh(sUm);
            uMesh.compute(ls);
            holeMesh = uMesh.createInnerMesh();
            
        end
  
        function ls = computeCircleLevelSet(obj, mesh,radius)
            gPar.type          = 'Circle';
            gPar.radius        = radius;
            gPar.xCoorCenter   = 0.5;
            gPar.yCoorCenter   = 0.5;
            g                  = GeometricalFunction(gPar);
            phiFun             = g.computeLevelSetFunction(mesh);
            lsCircle           = phiFun.fValues;
            ls = -lsCircle;
        end

        function s = updateCoordsMesh(obj, s)
            % Nudge nodes at the four rectangle corners in x to avoid
            % exact coincidences.
            tol  = 1e-8;
            epsx = 1e-9;
        
            x = s.coord(:,1); y = s.coord(:,2);
            xmax = max(x); xmin = min(x);
            ymax = max(y); ymin = min(y);
        
            % Top-right (xmax,ymax)
            mask = abs(x - xmax) < tol & abs(y - ymax) < tol;
            s.coord(mask, :) = s.coord(mask, :) - [epsx, 0];
        
            % Bottom-right (xmax,ymin)
            mask = abs(x - xmax) < tol & abs(y - ymin) < tol;
            s.coord(mask, :) = s.coord(mask, :) - [epsx, 0];
        
            % Top-left (xmin,ymax)
            mask = abs(x - xmin) < tol & abs(y - ymax) < tol;
            s.coord(mask, :) = s.coord(mask, :) + [epsx, 0];
        
            % Bottom-left (xmin,ymin)
            mask = abs(x - xmin) < tol & abs(y - ymin) < tol;
            s.coord(mask, :) = s.coord(mask, :) + [epsx, 0];
        end
      

        function [mD,mSb,iC,lG,iCR,discMesh] = createMeshDomain(obj)
            s.nsubdomains   = obj.nSubdomains;
            s.meshReference = obj.referenceMesh;
            s.tolSameNode = obj.tolSameNode;
            m = MeshCreatorFromRVE.create(s);
            [mD,mSb,iC,~,lG,iCR,discMesh] = m.create();
        end

        
        function mCoarse = createCoarseMesh(obj,mR)
            s.nsubdomains   = obj.nSubdomains; %nx ny
            s.meshReference = obj.createReferenceCoarseMesh(mR);
            s.tolSameNode   = obj.tolSameNode;
            mRVECoarse      = MeshCreatorFromRVE.create(s);
            [mCoarse,~,~] = mRVECoarse.create();
        end

        function cMesh = createReferenceCoarseMesh(obj,mR)
            xmax = max(mR.coord(:,1));
            xmin = min(mR.coord(:,1));
            ymax = max(mR.coord(:,2));
            ymin = min(mR.coord(:,2));
            coord(1,1) = xmin;  coord(1,2) = ymin;
            coord(2,1) = xmax;  coord(2,2) = ymin;
            coord(3,1) = xmax;  coord(3,2) = ymax;
            coord(4,1) = xmin;  coord(4,2) = ymax;
            connec = [2 3 4 1];
            s.coord = coord;
            s.connec = connec;
            cMesh = Mesh.create(s);
        end

        function createBCapplier(obj)
            s.mesh                  = obj.meshDomain;
            s.boundaryConditions    = obj.boundaryConditions;
            obj.bcApplier           = BCApplier(s);
        end

        function material = createMaterial(obj,mesh)
            E  = 1;
            nu = 1/3;  
            s.type    = 'ISOTROPIC';
            s.ptype   = 'ELASTIC';
            s.ndim    = mesh.ndim;
            s.young   = ConstantFunction.create(E,mesh);
            s.poisson = ConstantFunction.create(nu,mesh);
            tensor    = Material.create(s);
            material  = tensor;
        end


        function [bC,Dir] = createBoundaryConditions(obj)
            minx = min(obj.meshDomain.coord(:,1));
            maxx = max(obj.meshDomain.coord(:,1));
            miny = min(obj.meshDomain.coord(:,2));
            maxy = max(obj.meshDomain.coord(:,2));
            tolBound = obj.tolSameNode;
            isLeft   = @(coor) (abs(coor(:,1) - minx)   < tolBound);
            isRight  = @(coor) (abs(coor(:,1) - maxx)   < tolBound);
            isBottom = @(coor) (abs(coor(:,2) - miny)   < tolBound);
            isTop    = @(coor) (abs(coor(:,2) - maxy)   < tolBound);
%             Dir{1}.domain    = @(coor) isLeft(coor);%| isRight(coor) ;
%             Dir{1}.direction = [1,2];
%             Dir{1}.value     = 0;
%             dirichletFun = DirichletCondition(obj.meshDomain, Dir{1});
% 
            mesh = obj.meshDomain;
%             PL.domain    = @(coor) isRight(coor);
%             PL.direction = 2;
%             PL.value     = -0.1;
%             pointload = TractionLoad(mesh,PL,'DIRAC');

             Dir{1}.domain    = @(coor) isLeft(coor);%| isRight(coor) ;
            Dir{1}.direction = [1,2];
            Dir{1}.value     = 0;

                        Dir{2}.domain    = @(coor) isRight(coor) ;
                        Dir{2}.direction = [1,2];
                        Dir{2}.value     = 0;
            dirichletFun=[];        
            for i = 1:numel(Dir)
                dir = DirichletCondition(obj.meshDomain, Dir{i});
                dirichletFun = [dirichletFun, dir];
            end

            PL.domain    = @(coor) isTop(coor);
            PL.direction = [2];
            PL.value     = [-0.1];
            pointload = TractionLoad(obj.meshDomain,PL,'DIRAC');

%             mesh = obj.meshDomain;
%             PL.domain    = @(coor) isRight(coor);
%             PL.direction = 2;
%             PL.value     = -0.1;
%             pointload = PointLoad(mesh,PL);
%             % need this because force applied in the face not in a point
%             pointload.values        = pointload.values/size(pointload.dofs,1);
%             fvalues                 = zeros(mesh.nnodes*mesh.ndim,1);
%             fvalues(pointload.dofs) = pointload.values;
%             fvalues                 = reshape(fvalues,mesh.ndim,[])';
%             pointload.fun.setFValues(fvalues);

            s.pointloadFun = pointload;
            s.dirichletFun = dirichletFun;
            s.periodicFun  =[];
            s.mesh         = mesh;
            bC             = BoundaryConditions(s);                        
        end

  

        function [LHSr,Mr, RHSr] = createElasticProblem(obj)
            u = LagrangianFunction.create(obj.meshDomain,obj.meshDomain.ndim,'P1');
            material = obj.createMaterial(obj.meshDomain);
            [lhs,LHSr] = obj.computeStiffnessMatrix(obj.meshDomain,u,material);
            [M,Mr] = obj.computeMassMatrix(obj.meshDomain,u);
            RHSr       = obj.computeForces(lhs,u);
        end

        function [LHS,LHSr] = computeStiffnessMatrix(obj,mesh,dispFun,C)
            % s.type     = 'ElasticStiffnessMatrix';
            % s.mesh     = mesh;
            % s.test     = dispFun;
            % s.trial    = dispFun;
            % s.material = mat;
            % s.quadratureOrder = 2;
            % lhs = LHSIntegrator.create(s);
            % LHS = lhs.compute();

            LHS = IntegrateLHS(@(u,v) DDP(SymGrad(v),DDP(C,SymGrad(u))),dispFun,dispFun,mesh,'Domain',2);
            LHSr = obj.bcApplier.fullToReducedMatrixDirichlet(LHS);
        end

       function RHS =  computeForces(obj,stiffness,u)
            bc  = obj.boundaryConditions;
            t   = bc.tractionFun;
            rhs = zeros(u.nDofs,1);
            if ~isempty(t)
                for i = 1:numel(t)
                    rhsi = t(i).computeRHS(u);
                    rhs  = rhs + rhsi;
                end
            end
            if strcmp(obj.solverType,'REDUCED')
                bc      = obj.boundaryConditions;
                dirich  = bc.dirichlet_dofs;
                dirichV = bc.dirichlet_vals;
                if ~isempty(dirich)
                    R = -stiffness(:,dirich)*dirichV;
                else
                    R = zeros(sum(obj.uFun.nDofs(:)),1);
                end
                rhs = rhs+R;
            end
             RHS = obj.bcApplier.fullToReducedVectorDirichlet(rhs);
       end

        function [Meifem,Kcoarse, Mcoarse,eifem, ss] = createEIFEMPreconditioner(obj,dir,iC,lG,bS,iCR,dMesh,radiusMesh)
            mR = obj.referenceMesh;
%           
            s.RVE           = TrainedRVE(EIFEoper);
            s.mesh          = obj.createCoarseMesh(mR);
%            s.mesh          = obj.loadCoarseMesh(mR);
            s.DirCond       = dir;
            s.nSubdomains = obj.nSubdomains;
            s.ddDofManager = obj.createDomainDecompositionDofManager(iC,lG,bS,mR,iCR);
            eifem           = EIFEM(s);

            ss.ddDofManager = obj.createDomainDecompositionDofManager(iC,lG,bS,mR,iCR);
            ss.EIFEMsolver = eifem;
            ss.bcApplier = obj.bcApplier;
            ss.dMesh     = dMesh;
            ss.type = 'EIFEM';
            eP = Preconditioner.create(ss);
            Meifem = @(r) eP.apply(r);
            u = LagrangianFunction.create(s.mesh, s.mesh.ndim,'P1');

            Kcoarse = repmat(EIFEoper.Kcoarse,[1,1,s.mesh.nelem]);
            Kcoarse = eifem.assembleMatrix(Kcoarse,u,u);
            Kcoarse = eifem.reduceMatrix(Kcoarse);
            Mcoarse = repmat(EIFEoper.Mcoarse,[1,1,s.mesh.nelem]);
            Mcoarse = eifem.assembleMatrix(Mcoarse,u,u);
            Mcoarse = eifem.reduceMatrix(Mcoarse);
        end

        function offlineTrainer(obj, mesh)

            data = Training(mesh);
            p = OfflineDataProcessor(data); % i don't want to have to run this if I use NN
            EIFEoper = p.computeROMbasis(radiusMesh);

        end
        
        function d = createDomainDecompositionDofManager(obj,iC,lG,bS,mR,iCR)
            s.nSubdomains     = obj.nSubdomains;
            s.interfaceConnec = iC;
            s.interfaceConnecReshaped = iCR;
            s.locGlobConnec   = lG;
            s.nBoundaryNodes  = bS{1}.mesh.nnodes;
            s.nReferenceNodes = mR.nnodes;
            s.nNodes          = obj.meshDomain.nnodes;
            s.nDimf           = obj.meshDomain.ndim;
            d = DomainDecompositionDofManager(s);
        end

        function Milu = createILUpreconditioner(obj,LHS)
            s.LHS = LHS;
            s.type = 'ILU';
            M = Preconditioner.create(s);
            Milu = @(r) M.apply(r);
        end

        function plotResidual(obj,residualPCG,errPCG,errAnormPCG,residualCG,errCG,errAnormCG)
            figure
            plot(residualPCG,'linewidth',2)
            hold on
            plot(residualCG,'linewidth',2)
            set(gca, 'YScale', 'log')
            legend({'CG + ILU-EIFEM-ILU','CG'},'FontSize',12)
            xlabel('Iteration')
            ylabel('Residual')

            figure
            plot(errPCG,'linewidth',2)
            hold on
            plot(errCG,'linewidth',2)
            set(gca, 'YScale', 'log')
            legend('CG + EIFEM+ ILU(CG-90%-L2)','CG')
            xlabel('Iteration')
            ylabel('||error||_{L2}')

            figure
            plot(errAnormPCG,'linewidth',2)
            hold on
            plot(errAnormCG,'linewidth',2)
            set(gca, 'YScale', 'log')
            legend('CG + EIFEM+ ILU(CG-90%-L2)','CG')
            xlabel('Iteration')
            ylabel('Energy norm')

        end
       
        function printContinuousSolution(obj,EIFEM,phiCoarse,mesh)
            for i=1:size(phiCoarse,2)
                PhiFineDisc(:,:,i) = EIFEM.reconstructSolution(phiCoarse(:,i));
                PhiFineCont(:,:,i) = EIFEM.computeContinousField(PhiFineDisc(:,:,i));  
                
                uplot = PhiFineCont(:,i);
                EIFEMtesting.plotSolution(uplot,mesh,15,2,i,[],0)
                %PhiFineContRed(:,:,i) = obj.bcApplier.fullToReducedVectorDirichlet(PhiFineCont(:,:,i));
                %save('PhiCoarseProjectedRed.mat','PhiFineContRed')

            end
        end

        function printDiscontinuousSolution(obj,EIFEM, phiCoarse,mesh)
            for i=1:size(phiCoarse,2)
                PhiFineDisc(:,:,i) = EIFEM.reconstructSolution(phiCoarse(:,i));
                uplot = PhiFineDisc(:,:,i);
                uplot = uplot(:);
                EIFEMtesting.plotSolution(uplot,mesh,15,2,i,[],0)
            end
        end

        function printCoarseSolution(obj,EIFEM, phiCoarse,mesh)
            for i=1:size(phiCoarse,2)
                PhiBoundary(:,:,i) = EIFEM.injectCoarseToFineBoundary(phiCoarse(:,i));
                uplot = PhiFineCont(:,i);
                %uplot = uplot(:);
                EIFEMtesting.plotSolution(uplot,mesh,15,2,i,[],0)
            end
        end

        function printFineSolution(obj,phiFine,mesh)
             for i=1:size(phiFine,2)
                 PhiFineFull(:,:,i) = obj.bcApplier.reducedToFullVectorDirichlet(phiFine(:,i));
                 uplot = PhiFineFull(:,:,i);
                 uplot = uplot(:);
                 EIFEMtesting.plotSolution(uplot,mesh,15,2,i,[],0)
            end
        end
    end

end
