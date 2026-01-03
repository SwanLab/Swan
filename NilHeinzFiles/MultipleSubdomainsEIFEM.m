classdef MultipleSubdomainsEIFEM < handle

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
        subdomainMeshes
        r
        xmin
        xmax
        ymin
        ymax
        Nr
        Ntheta
        x0
        y0
        iC
    end


    properties (Access = private)
        nSubdomains
    end

    methods (Access = public)

        function obj = MultipleSubdomainsEIFEM()
            close all
            obj.init()

            [mD,mSb,iC,lG,iCR,discMesh,bS] = obj.createMesh();
            obj.meshDomain = mD;
            obj.subdomainMeshes = mSb;
            

            [bC,dir] = obj.createBoundaryConditions();
            obj.boundaryConditions = bC;
            obj.createBCapplier()

            [LHSr, Mr, RHSr] = obj.createElasticProblem();

            LHSfun = @(x) LHSr*x;
            
            [Meifem, EIFEM, ss]       = obj.createEIFEMPreconditioner(dir,iC,lG,bS,iCR,discMesh,obj.r);

            Milu         = obj.createILUpreconditioner(LHSr);
            Mmult        = @(r) Preconditioner.multiplePrec(r,LHSfun,Milu,Meifem,Milu);
            Mid          = @(r) r;

            tol = 1e-8;
            x0 = zeros(size(RHSr));
            xSol = LHSr\RHSr;

            [uPCG,residualPCG,errPCG,errAnormPCG] = PCG.solve(LHSfun,RHSr,x0,Mmult,tol,xSol);     
            [uCG,residualCG,errCG,errAnormCG]    = PCG.solve(LHSfun,RHSr,x0,Mid,tol,xSol);

            obj.plotResidual(residualPCG,errPCG,errAnormPCG,residualCG,errCG,errAnormCG)

            uCGFull = obj.bcApplier.reducedToFullVectorDirichlet(uCG); %reapply all dirichlett DOFs
            uF = saveDeformed(obj.meshDomain,uCGFull);
            plot(uF)

        end

    end

    methods (Access = private)

        function init(obj)
            obj.nSubdomains = [5,3];
            obj.r = 1e-6*ones(obj.nSubdomains)'; 
            obj.r= (0.6 - 0.2) * rand(obj.nSubdomains(2),obj.nSubdomains(1)) + 0.1;
            obj.xmax=1; obj.xmin=-1; obj.ymax = 1; obj.ymin=-1;
            obj.Nr = 7; obj.Ntheta = 14; 
            obj.x0 = 0; obj.y0=0;
            obj.tolSameNode = 1e-10;
            obj.solverType = 'REDUCED';
        end  

        function [mD,mSb,iC,lG,iCR,discMesh,bS] = createMesh(obj)
            mSbd = obj.createSubDomainMeshes();
            bS = mSbd{1,1}.createBoundaryMesh();
            [mD,mSb,iC,lG,iCR,discMesh] = obj.createMeshDomainJoiner(mSbd);
            
        end

         function  mSbd = createSubDomainMeshes(obj)
            nX = obj.nSubdomains(1);
            nY = obj.nSubdomains(2);
            Lx = obj.xmax-obj.xmin;
            Ly = obj.ymax-obj.ymin;
            for jDom = 1:nY
                for iDom = 1:nX
                    refMesh = mesh_rectangle_via_triangles(obj.r(jDom,iDom),obj.xmax,obj.xmin,obj.ymax,obj.ymin,obj.Nr,obj.Ntheta,obj.x0,obj.y0);
                    coord0 = refMesh.coord;
                    s.coord(:,1) = coord0(:,1)+Lx*(iDom-1);
                    s.coord(:,2) = coord0(:,2)+Ly*(jDom-1);
                    s.connec = refMesh.connec;
                    mIJ     = Mesh.create(s);
                    %                     plot(mIJ)
                    %                     hold on;
                    mSbd{jDom,iDom} = mIJ;
                end
            end
            obj.referenceMesh = mSbd{1,1};
        end

         function [mD,mSb,iC,lG,iCR,discMesh] = createMeshDomainJoiner(obj,mSbd)
            s.nsubdomains   = obj.nSubdomains;
            s.meshReference = obj.referenceMesh;
            s.tolSameNode = obj.tolSameNode;
            s.meshSbd     = mSbd;
            %             m = MeshCreatorFromRVE.create(s);
            m = MeshJoiner(s);
            [mD,mSb,iC,~,lG,iCR,discMesh] = m.create();
         end

        function createReferenceMesh(obj)

            squareMesh  = UnitTriangleMesh(12,12);
            s.coord     = squareMesh.coord;
            s.connec    = squareMesh.connec;
            s.interType = 'LINEAR';
            s           = obj.updateCoordsMesh(s);
            obj.referenceMesh = Mesh.create(s);
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
            tolBound = obj.tolSameNode;
            isLeft   = @(coor) (abs(coor(:,1) - minx)   < tolBound);
            isRight  = @(coor) (abs(coor(:,1) - maxx)   < tolBound);
            Dir{1}.domain    = @(coor) isLeft(coor);%| isRight(coor) ;
            Dir{1}.direction = [1,2];
            Dir{1}.value     = 0;
            dirichletFun = DirichletCondition(obj.meshDomain, Dir{1});

            mesh = obj.meshDomain;
            PL.domain    = @(coor) isRight(coor);
            PL.direction = 2;
            PL.value     = -0.1;
            pointload = TractionLoad(mesh,PL,'DIRAC');
            % need this because force applied in the face not in a point
            % pointload.values        = pointload.values/size(pointload.dofs,1);
            % fvalues                 = zeros(mesh.nnodes*mesh.ndim,1);
            % fvalues(pointload.dofs) = pointload.values;
            % fvalues                 = reshape(fvalues,mesh.ndim,[])';
            % pointload.fun.setFValues(fvalues);

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

        function [M , Mr] = computeMassMatrix(obj, mesh, dispFun)
            rho = obj.computeDensity(mesh);
            M = IntegrateLHS(@(u,v) rho .* DP(v,u),dispFun,dispFun,mesh,'Domain',2); %DP(u,v) 
            Mr = obj.bcApplier.fullToReducedMatrixDirichlet(M);

        end

        function rho = computeDensity(obj, mesh)
            rho = 1; % kg/m^3 -
            rho  = ConstantFunction.create(rho, mesh);
        end

        function [lambda, Phi, omega] = computeModalAnalysis(obj, K, M)
            [Phi, D] = eigs(K, M, 15, "smallestabs");
            lambda = diag(D);
            [lambda, idx] = sort(lambda, 'ascend'); 
            Phi = Phi(:, idx); %sort
            omega = sqrt(max(lambda,0));
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

        function [Meifem,eifem, ss] = createEIFEMPreconditioner(obj,dir,iC,lG,bS,iCR,dMesh,radiusMesh)
            mR = obj.referenceMesh;
  
            s.RVE           = TrainedRVEMultipleRadius(radiusMesh,mR);
            s.mesh          = obj.createCoarseMesh(mR);
            s.DirCond       = dir;
            s.nSubdomains = obj.nSubdomains;
            s.ddDofManager = obj.createDomainDecompositionDofManager(iC,lG,bS,mR,iCR);
            eifem           = EIFEMdifferentSubd(s);

            ss.ddDofManager = obj.createDomainDecompositionDofManager(iC,lG,bS,mR,iCR);
            ss.EIFEMsolver = eifem;
            ss.bcApplier = obj.bcApplier;
            ss.dMesh     = dMesh;
            ss.type = 'EIFEM';
            eP = Preconditioner.create(ss);
            Meifem = @(r) eP.apply(r);

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
