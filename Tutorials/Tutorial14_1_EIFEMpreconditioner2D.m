classdef Tutorial14_1_EIFEMpreconditioner2D < handle

    properties (Access = private)
        referenceMesh
        mesh
        boundaryConditions
        bcApplier
        LHS
        RHS
        fileNameEIFEM
        tolSameNode
        solverType
        young
        poisson
        EIFEMdata
        tolResidual
        nSubdomains
    end

    methods (Access = public)

        function obj = Tutorial14_1_EIFEMpreconditioner2D()
            close all
            obj.init()
            
            [obj.mesh,mSb,iC,lG,iCR,discMesh,bS] = obj.createMesh();
            [obj.boundaryConditions,dir]         = obj.createBoundaryConditions();
            obj.bcApplier                        = obj.createBCapplier();
            [LHSr,RHSr]                          = obj.createElasticProblem();
            Pmult                                = obj.createPreconditioner(dir,iC,lG,bS,iCR,discMesh,LHSr);
            pcg                                  = obj.createSolver(Pmult,zeros(size(RHSr)),obj.tolResidual);
            [uPCG,residualPCG,~,~]               = pcg.solve(LHSr,RHSr);     
           
            % xSol = LHSr\RHSr;
            % [uPCG,residualPCG,errPCG,errAnormPCG] = pcg.solve(LHSr,RHSr,xSol);     
            % [uCG,residualCG,errCG,errAnormCG]    = PCG.solve(LHSfun,RHSr,x0,Mid,tol,xSol);     
            % obj.plotResidual(residualPCG,errPCG,errAnormPCG,residualCG,errCG,errAnormCG)
        end

    end

    methods (Access = private)

        function init(obj)
            obj.nSubdomains             = [15 5]; %nx ny
            obj.fileNameEIFEM           = 'DEF_Q4porL_1.mat';
            obj.tolSameNode             = 1e-10;
            obj.solverType              = 'REDUCED';
            obj.EIFEMdata.nSubdomains   = [5,5];
            obj.EIFEMdata.domainIndices = [3,3];
            obj.EIFEMdata.tolSameNode   = 1e-11;
            obj.EIFEMdata.Coarseorder   = 1;
            obj.EIFEMdata.type          = 'Continuous';
            obj.tolResidual             = 1e-8;
        end        

        function createReferenceMesh(obj)
            filename = obj.fileNameEIFEM;
            load(filename);
            s.coord    = EIFEoper.MESH.COOR;
            s.connec   = EIFEoper.MESH.CN;
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

        function [mD,mSb,iC,lG,iCR,discMesh,bS] = createMesh(obj)
             obj.createReferenceMesh();
             bS = obj.referenceMesh.createBoundaryMesh();
             [mD,mSb,iC,lG,iCR,discMesh] = obj.createMeshDomain();
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
            connec   = [1 2 3 4];
            s.coord  = coord;
            s.connec = connec;
            cMesh = Mesh.create(s);
        end

        function bcApplier = createBCapplier(obj)
            s.mesh               = obj.mesh;
            s.boundaryConditions = obj.boundaryConditions;
            bcApplier            = BCApplier(s);
        end

        function material = createMaterial(obj,mesh)
            obj.young = 1;
            obj.poisson    = 1/3;  

            N = obj.mesh.ndim;
            E = ConstantFunction.create(obj.young,mesh);
            nu = ConstantFunction.create(obj.poisson,mesh); 
            lam = LameParametersConverter;

            mu     = lam.computeShearFromYoungAndPoisson(E,nu);  mu = Expand(mu,4);
            lambda = lam.computeLambdaFromYoungAndPoisson(E,nu,N);  lambda = Expand(lambda,4);
            I      = ConstantFunction.create(eye4D(N),obj.mesh);
            IxI    = ConstantFunction.create(kronEye(N),obj.mesh);
            
            material = 2*mu.*I + lambda.*IxI;    
        end

        function [bC,Dir] = createBoundaryConditions(obj)
            minx = min(obj.mesh.coord(:,1));  maxx = max(obj.mesh.coord(:,1));
            miny = min(obj.mesh.coord(:,2));  maxy = max(obj.mesh.coord(:,2));
            tolBound = obj.tolSameNode;
            isLeft   = @(coor) (abs(coor(:,1) - minx)   < tolBound);
            isRight  = @(coor) (abs(coor(:,1) - maxx)   < tolBound);
            isBottom = @(coor) (abs(coor(:,2) - miny)   < tolBound);
            isTop    = @(coor) (abs(coor(:,2) - maxy)   < tolBound);
 
            mesh = obj.mesh;

            Dir{1}.domain    = @(coor) isLeft(coor);%| isRight(coor) ;
            Dir{1}.direction = [1,2];
            Dir{1}.value     = 0;

            Dir{2}.domain    = @(coor) isRight(coor) ;
            Dir{2}.direction = [1,2];
            Dir{2}.value     = 0;
            dirichletFun=[]; 

            for i = 1:numel(Dir)
                dir = DirichletCondition(obj.mesh, Dir{i});
                dirichletFun = [dirichletFun, dir];
            end

            PL.domain    = @(coor) isTop(coor);
            PL.direction = [2];
            PL.value     = [-0.1];
            pointload = TractionLoad(obj.mesh,PL,'DIRAC');


            s.pointloadFun = pointload;
            s.dirichletFun = dirichletFun;
            s.periodicFun  =[];
            s.mesh         = mesh;
            bC             = BoundaryConditions(s);                        
        end  

        function [LHSr,RHSr] = createElasticProblem(obj)
            u          = LagrangianFunction.create(obj.mesh,obj.mesh.ndim,'P1');
            material   = obj.createMaterial(obj.mesh);
            [lhs,LHSr] = obj.computeStiffnessMatrix(obj.mesh,u,material);
            RHSr       = obj.computeForces(lhs,u);
        end

        function [LHS,LHSr] = computeStiffnessMatrix(obj,mesh,dispFun,C)
            LHS  = IntegrateLHS(@(u,v) DDP(SymGrad(v),DDP(C,SymGrad(u))),dispFun,dispFun,mesh,'Domain',2);
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
            ss.ddDofManager = obj.createDomainDecompositionDofManager(iC,lG,bS,obj.referenceMesh,iCR);
            ss.EIFEMsolver  = obj.createEIFEM(dir);
            ss.bcApplier    = obj.bcApplier;
            ss.dMesh        = dMesh;
            ss.type         = 'EIFEM';
            Meifem = Preconditioner.create(ss);
            % Meifem = @(r) eP.apply(r);
        end

        function EIFEoper = computeEIFEMoperators(obj)
            data     = obj.computeTrainingData();
            EIFEoper = obj.computeReducedOperators(data);         
        end

        function data = computeTrainingData(obj)
            s         = obj.EIFEMdata;
            s.mesh    = obj.referenceMesh;
            s.young   = obj.young;
            s.poisson = obj.poisson;
            training = EIFEMTraining(s);
            data     = training.train();
        end

        function EIFEoper = computeReducedOperators(obj,data)
            data.Coarseorder = obj.EIFEMdata.Coarseorder;
            data.type        = obj.EIFEMdata.type;
            data.young       = obj.young;
            data.poisson     = obj.poisson;
            p        = OfflineDataProcessor(data);
            EIFEoper = p.computeROMbasis();
        end
        
        function eifem = createEIFEM(obj,dir)
            EIFEoper      = obj.computeEIFEMoperators();
            s.RVE         = TrainedRVE(EIFEoper);
            s.mesh        = obj.createCoarseMesh(obj.referenceMesh);
%            s.mesh          = obj.loadCoarseMesh(mR);
            s.DirCond     = dir;
            s.nSubdomains = obj.nSubdomains;
            eifem         = EIFEM(s);
        end
        
        function d = createDomainDecompositionDofManager(obj,iC,lG,bS,mR,iCR)
            s.nSubdomains             = obj.nSubdomains;
            s.interfaceConnec         = iC;
            s.interfaceConnecReshaped = iCR;
            s.locGlobConnec           = lG;
            s.nBoundaryNodes          = bS{1}.mesh.nnodes;
            s.nReferenceNodes         = mR.nnodes;
            s.nNodes                  = obj.mesh.nnodes;
            s.nDimf                   = obj.mesh.ndim;
            d = DomainDecompositionDofManager(s);
        end

        function Milu = createILUpreconditioner(obj,LHS)
            s.LHS  = LHS;
            s.type = 'ILU';
            Milu = Preconditioner.create(s);
        end

        function M = createMultiplicativePreconditioner(obj,varargin)
            s.varargin = varargin;
            s.type     = 'Multiplicative';
            M = Preconditioner.create(s);
        end

        function Mmult = createPreconditioner(obj,dir,iC,lG,bS,iCR,discMesh,LHSr)
            Meifem = obj.createEIFEMPreconditioner(dir,iC,lG,bS,iCR,discMesh);
            Milu   = obj.createILUpreconditioner(LHSr);
            Mmult  = obj.createMultiplicativePreconditioner(Milu,Meifem,Milu);
        end

        function solver = createSolver(obj,P,x0,tol)
            s.x0             = x0;            
            s.preconditioner = P;
            s.tol            = tol;
            solver = PCG(s);
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
