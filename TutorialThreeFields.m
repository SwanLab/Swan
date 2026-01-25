classdef TutorialThreeFields < handle

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

        function obj = TutorialThreeFields()
            close all
            obj.init()

            obj.createReferenceMesh();
            bS  = obj.referenceMesh.createBoundaryMesh();
            [mD,mSb,iC,lG,iCR,discMesh] = obj.createMeshDomain();
            obj.meshDomain = mD;
            mD.plot()
            [bC,dir] = obj.createBoundaryConditions();
            obj.boundaryConditions = bC;
            obj.createBCapplier()

            [LHSr,RHSr] = obj.createElasticProblem();

            LHSfun = @(x) LHSr*x;
            Meifem       = obj.createEIFEMPreconditioner(dir,iC,lG,bS,iCR,discMesh);
            Milu         = obj.createILUpreconditioner(LHSr);
            Mmult        = @(r) Preconditioner.multiplePrec(r,LHSfun,Milu,Meifem,Milu);
            Mid          = @(r) r;


            tol = 1e-8;
            x0 = zeros(size(RHSr));
            xSol = LHSr\RHSr;

            [uPCG,residualPCG,errPCG,errAnormPCG] = PCG.solve(LHSfun,RHSr,x0,Mmult,tol,xSol,obj.meshDomain,obj.bcApplier);     
            [uCG,residualCG,errCG,errAnormCG]    = PCG.solve(LHSfun,RHSr,x0,Mid,tol,xSol);     
            obj.plotResidual(residualPCG,errPCG,errAnormPCG,residualCG,errCG,errAnormCG)
        end

    end

    methods (Access = private)

        function init(obj)
            obj.nSubdomains  = [15 5]; %nx ny
            obj.fileNameEIFEM = 'DEF_Q4porL_1.mat';
            obj.tolSameNode = 1e-10;
            obj.solverType = 'REDUCED';
        end        

        function createReferenceMesh(obj)
            filename = obj.fileNameEIFEM;
            load(filename);
            s.coord    = EIFEoper.MESH.COOR;
            s.connec   = EIFEoper.MESH.CN;
            s.interType = 'QUADRATIC';
            obj.referenceMesh = Mesh.create(s);

            %% Uncomment for meshes that have corners and generate the mesh with the updated coordinates
%             tol = 1e-8;
%             xmax = max(s.coord(:,1)); xmin = min(s.coord(:,1));
%             ymax = max(s.coord(:,2)); ymin = mmin(s.coord(:,2));
%              % Top-right corner (xmax, ymax)
%             mask = abs(s.coord(:,1) - xmax) < tol & abs(s.coord(:,2) - ymax) < tol;
%             s.coord(mask, :) = s.coord(mask, :) - [1e-9, 0];
% 
%             % Bottom-right corner (xmax, ymin)
%             mask = abs(s.coord(:,1) - xmax) < tol & abs(s.coord(:,2) - ymin) < tol;
%             s.coord(mask, :) = s.coord(mask, :) - [1e-9, 0];
% 
%             % Top-left corner (xmin, ymax)
%             mask = abs(s.coord(:,1) - xmin) < tol & abs(s.coord(:,2) - ymax) < tol;
%             s.coord(mask, :) = s.coord(mask, :) + [1e-9, 0];
% 
%             % Bottom-left corner (xmin, ymin)
%             mask = abs(s.coord(:,1) - xmin) < tol & abs(s.coord(:,2) - ymin) < tol;
%             s.coord(mask, :) = s.coord(mask, :) + [1e-9, 0];
%             obj.referenceMesh = Mesh.create(s);
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

  

        function [LHSr,RHSr] = createElasticProblem(obj)
            u = LagrangianFunction.create(obj.meshDomain,obj.meshDomain.ndim,'P1');
            material = obj.createMaterial(obj.meshDomain);
            [lhs,LHSr] = obj.computeStiffnessMatrix(obj.meshDomain,u,material);
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

        function Meifem = createEIFEMPreconditioner(obj,dir,iC,lG,bS,iCR,dMesh)
            mR = obj.referenceMesh;
%             % obj.EIFEMfilename = '/home/raul/Documents/Thesis/EIFEM/RAUL_rve_10_may_2024/EXAMPLE/EIFE_LIBRARY/DEF_Q4porL_2s_1.mat';
%             EIFEMfilename = obj.fileNameEIFEM;
%             % obj.EIFEMfilename = '/home/raul/Documents/Thesis/EIFEM/05_HEXAG2D/EIFE_LIBRARY/DEF_Q4auxL_1.mat';
%             filename        = EIFEMfilename;
%             s.RVE           = TrainedRVE(filename);
            data = Training(mR);
            p = OfflineDataProcessor(data);
            EIFEoper = p.computeROMbasis();
            s.RVE           = TrainedRVE(EIFEoper);
            s.mesh          = obj.createCoarseMesh(mR);
%            s.mesh          = obj.loadCoarseMesh(mR);
            s.DirCond       = dir;
            s.nSubdomains = obj.nSubdomains;
            eifem           = EIFEM(s);

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
       
    end

end
