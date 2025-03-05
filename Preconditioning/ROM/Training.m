classdef Training < handle

    properties (Access = public)

    end
    properties (Access = private)
        meshDomain
        boundaryMeshJoined
        localGlobalConnecBd
        DirFun
        boundaryConditions
        bcApplier
        LHS
        RHS

        fileNameEIFEM
        tolSameNode

    end


    properties (Access = private)
        nSubdomains
    end

    methods (Access = public)

        function obj = Training()
            close all
            obj.init()

            mR = obj.createReferenceMesh();
            bS  = mR.createBoundaryMesh();
            [mD,mSb,iC,lG,iCR,discMesh] = obj.createMeshDomain(mR);
            obj.meshDomain = mD;
            [obj.boundaryMeshJoined, obj.localGlobalConnecBd] = obj.meshDomain.createSingleBoundaryMesh();
            obj.DirFun = obj.AnalyticalDirCond();
%             [bC,dir] = obj.createBoundaryConditions(obj.meshDomain);
%             obj.boundaryConditions = bC;
%             obj.createBCapplier()

            [LHS,RHS,uFun,lambdaFun] = obj.createElasticProblem();
            sol = LHS\RHS;
            u   = sol(1:uFun.nDofs,:);
%             obj.LHS = LHSf;
%             %             LHS = 0.5*(LHS+LHS');
% 
%             LHSf = @(x) LHS*x;
%             RHSf = RHS;
%             Usol = LHS\RHS;
%             Ufull = obj.bcApplier.reducedToFullVectorDirichlet(Usol);
%             %obj.plotSolution(Ufull,obj.meshDomain,1,1,0,obj.bcApplier,0)
% 
% 
%             Mid          = @(r) r;
%             Meifem       = obj.createEIFEMPreconditioner(mR,dir,iC,lG,bS,iCR,discMesh);
%             Milu         = obj.createILUpreconditioner(LHS);
%             MgaussSeidel = obj.createGaussSeidelpreconditioner(LHS);
%             MJacobi      = obj.createJacobipreconditioner(LHS);
%             Mmodal       = obj.createModalpreconditioner(LHS);
%             %             MdirNeu      = obj.createDirichletNeumannPreconditioner(mR,dir,iC,lG,bS,obj.LHS,mSb);
% 
%             MiluCG = @(r,iter) Preconditioner.InexactCG(r,LHSf,Milu,RHSf);
% 
%             tol = 1e-8;
%             tic
%             x0 = zeros(size(RHSf));
% %             [uCG,residualCG,errCG,errAnormCG] = PCG.solve(LHSf,RHSf,x0,Milu,tol,Usol,obj.meshDomain,obj.bcApplier);
%             toc
%             %             [uCG,residualCG,errCG,errAnormCG] = RichardsonSolver.solve(LHSf,RHSf,x0,P,tol,0.1,Usol);
% 
%             tol = 1e-8;
% 
%             %Mmult = MdirNeu;
%             x0 = zeros(size(RHSf));
%             r = RHSf - LHSf(x0);
%             Mmult = @(r,uk) Preconditioner.multiplePrec(r,MiluCG,Meifem,MiluCG,LHSf,RHSf,obj.meshDomain,obj.bcApplier,uk);
% %              Mmult = @(r) Preconditioner.multiplePrec(r,Mid,Meifem,Mid,LHSf,RHSf,obj.meshDomain,obj.bcApplier);
% %             zmult = Mmult(r);
%             
% %             zfull = obj.bcApplier.reducedToFullVectorDirichlet(zmult);
%             %obj.plotSolution(zfull,obj.meshDomain,0,0,2,obj.bcApplier,0)
% 
% %             zeifem = Meifem(r);
% %             zfull = obj.bcApplier.reducedToFullVectorDirichlet(zeifem);
%             %obj.plotSolution(zfull,obj.meshDomain,0,0,1,obj.bcApplier,0)
%            % x0 = zmult;
%             tic
%             %           tau = @(r,A) 1;
%             [uPCG,residualPCG,errPCG,errAnormPCG] = PCG.solve(LHSf,RHSf,x0,Mmult,tol,Usol,obj.meshDomain,obj.bcApplier);
%             %            [uCG,residualPCG,errPCG,errAnormPCG] = RichardsonSolver.solve(LHSf,RHSf,x0,Mmult,tol,tau,Usol);
%             toc
% 
%             figure
%             plot(residualPCG,'linewidth',2)
%             hold on
%             plot(residualCG,'linewidth',2)
%             set(gca, 'YScale', 'log')
%             legend({'CG + ILU-EIFEM-ILU','CG'},'FontSize',12)
%             xlabel('Iteration')
%             ylabel('Residual')
% 
%             figure
%             plot(errPCG,'linewidth',2)
%             hold on
%             plot(errCG,'linewidth',2)
%             set(gca, 'YScale', 'log')
%             legend('CG + EIFEM+ ILU(CG-90%-L2)','CG')
%             xlabel('Iteration')
%             ylabel('||error||_{L2}')
% 
%             figure
%             plot(errAnormPCG,'linewidth',2)
%             hold on
%             plot(errAnormCG,'linewidth',2)
%             set(gca, 'YScale', 'log')
%             legend('CG + EIFEM+ ILU(CG-90%-L2)','CG')
%             xlabel('Iteration')
%             ylabel('Energy norm')

        end

    end

    methods (Access = private)

        function init(obj)
            obj.nSubdomains  = [5 5]; %nx ny
%             obj.fileNameEIFEM = 'DEF_Q4auxL_1.mat';
            obj.fileNameEIFEM = 'DEF_Q4porL_1.mat';
            obj.tolSameNode = 1e-10;
        end

        function [mD,mSb,iC,lG,iCR,discMesh] = createMeshDomain(obj,mR)
            s.nsubdomains   = obj.nSubdomains; %nx ny
            s.meshReference = mR;
            s.tolSameNode = obj.tolSameNode;
            m = MeshCreatorFromRVE(s);
            [mD,mSb,iC,~,lG,iCR,discMesh] = m.create();
        end


        function mS = createReferenceMesh(obj)
            %     mS = obj.createStructuredMesh();
            %   mS = obj.createMeshFromGid();
            mS = obj.createEIFEMreferenceMesh();
        end


        function mS = createMeshFromGid(obj)
            filename   = 'lattice_ex1';
            a.fileName = filename;
            femD       = FemDataContainer(a);
            mS         = femD.mesh;
        end

        function mS = createStructuredMesh(obj)

            % Generate coordinates
            x1 = linspace(0,1,2);
            x2 = linspace(0,1,2);
            % Create the grid
            [xv,yv] = meshgrid(x1,x2);
            % Triangulate the mesh to obtain coordinates and connectivities
            [F,coord] = mesh2tri(xv,yv,zeros(size(xv)),'x');

            s.coord    = coord(:,1:2);
            s.connec   = F;
            mS         = Mesh.create(s);
        end

        function mS = createEIFEMreferenceMesh(obj)
            filename = obj.fileNameEIFEM;
            load(filename);
            s.coord    = EIFEoper.MESH.COOR;
            isMin = s.coord==min(s.coord);
            isMax = s.coord==max(s.coord);

            s.connec   = EIFEoper.MESH.CN;
            mS         = Mesh.create(s);
        end

        function mCoarse = createCoarseMesh(obj,mR)
            s.nsubdomains   = obj.nSubdomains; %nx ny
            s.meshReference = obj.createReferenceCoarseMesh(mR);
            s.tolSameNode   = obj.tolSameNode;
            mRVECoarse      = MeshCreatorFromRVE(s);
            [mCoarse,~,~] = mRVECoarse.create();
        end



        function cMesh = createReferenceCoarseMesh(obj,mR)
            xmax = max(mR.coord(:,1));
            xmin = min(mR.coord(:,1));
            ymax = max(mR.coord(:,2));
            ymin = min(mR.coord(:,2));
            coord(1,1) = xmin;
            coord(1,2) = ymin;
            coord(2,1) = xmax;
            coord(2,2) = ymin;
            coord(3,1) = xmax;
            coord(3,2) = ymax;
            coord(4,1) = xmin;
            coord(4,2) = ymax;
            %             coord(1,1) = xmax;
            %             coord(1,2) = ymin;
            %             coord(2,1) = xmax;
            %             coord(2,2) = ymax;
            %             coord(3,1) = xmin;
            %             coord(3,2) = ymax;
            %             coord(4,1) = xmin;
            %             coord(4,2) = ymin;
            connec = [1 2 3 4];
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
            [young,poisson] = obj.computeElasticProperties(mesh);
            s.type    = 'ISOTROPIC';
            s.ptype   = 'ELASTIC';
            s.ndim    = mesh.ndim;
            s.young   = young;
            s.poisson = poisson;
            tensor    = Material.create(s);
            material  = tensor;
        end

        function [young,poisson] = computeElasticProperties(obj,mesh)
            E  = 1;
            nu = 1/3;
            young   = ConstantFunction.create(E,mesh);
            poisson = ConstantFunction.create(nu,mesh);
        end

        function [Dir,PL] = createRawBoundaryConditions(obj)
            minx = min(obj.meshDomain.coord(:,1));
            maxx = max(obj.meshDomain.coord(:,1));
            miny = min(obj.meshDomain.coord(:,2));
            maxy = max(obj.meshDomain.coord(:,2));
            tolBound = obj.tolSameNode;
            isLeft   = @(coor) (abs(coor(:,1) - minx)   < tolBound);
            isRight  = @(coor) (abs(coor(:,1) - maxx)   < tolBound);
            isBottom = @(coor) (abs(coor(:,2) - miny)   < tolBound);
            isTop    = @(coor) (abs(coor(:,2) - maxy)   < tolBound);
            %             isMiddle = @(coor) (abs(coor(:,2) - max(coor(:,2)/2)) == 0);
            Dir{1}.domain    = @(coor) isLeft(coor);%| isRight(coor) ;
            Dir{1}.direction = [1,2];
            Dir{1}.value     = 0;

            %             Dir{2}.domain    = @(coor) isRight(coor) ;
            %             Dir{2}.direction = [2];
            %             Dir{2}.value     = 0;

%             PL.domain    = @(coor) isTop(coor);
%             PL.direction = [2];
%             PL.value     = [-0.1];
                        PL.domain    = @(coor) isRight(coor);
                        PL.direction = [1];
                        PL.value     = [0.1];
        end

        function [bc,Dir,PL] = createBoundaryConditions(obj,mesh)
            [Dir,PL]  = obj.createRawBoundaryConditions();
            dirichletFun = [];
            for i = 1:numel(Dir)
                dir = DirichletCondition(obj.meshDomain, Dir{i});
                dirichletFun = [dirichletFun, dir];
            end

            pointload = PointLoad(mesh,PL);
            % need this because force applied in the face not in a point
            pointload.values        = pointload.values/size(pointload.dofs,1);
            fvalues                 = zeros(mesh.nnodes*mesh.ndim,1);
            fvalues(pointload.dofs) = pointload.values;
            fvalues                 = reshape(fvalues,mesh.ndim,[])';
            pointload.fun.setFValues(fvalues);

            s.pointloadFun = pointload;
            s.dirichletFun = dirichletFun;
            s.periodicFun  =[];
            s.mesh         = mesh;
            bc             = BoundaryConditions(s);
        end


        function [LHS,RHS,u,dLambda] = createElasticProblem(obj)
            u = LagrangianFunction.create(obj.meshDomain,obj.meshDomain.ndim,'P1');
            dLambda = LagrangianFunction.create(obj.boundaryMeshJoined,obj.meshDomain.ndim,'P1');
            LHS = obj.computeLHS(u,dLambda);
            RHS = obj.computeRHS(u,dLambda);
%             RHS = obj.computeForces(lhs,u);
        end

        function LHS  = computeLHS(obj,u,dLambda)
           
            material = obj.createMaterial(obj.meshDomain);
            K = obj.computeStiffnessMatrix(obj.meshDomain,u,material);
           
            C = obj.computeConditionMatrix(dLambda);
            Z = zeros(size(C,2));
            LHS = [K C; C' Z];
        end

        function LHS = computeStiffnessMatrix(obj,mesh,dispFun,mat)
            s.type     = 'ElasticStiffnessMatrix';
            s.mesh     = mesh;
            s.test     = dispFun;
            s.trial    = dispFun;
            s.material = mat;
            s.quadratureOrder = 2;
            lhs = LHSIntegrator.create(s);
            LHS = lhs.compute();
%             LHSr = obj.bcApplier.fullToReducedMatrixDirichlet(LHS);
        end

        function C = computeConditionMatrix(obj,dLambda)
            s.type                  = 'ConditionMatrix';   
            s.quadType              = 2;
            s.boundaryMeshJoined    = obj.boundaryMeshJoined;
            s.localGlobalConnecBd   = obj.localGlobalConnecBd;
            s.nnodes                = obj.meshDomain.nnodes;
            lhs      = LHSIntegrator_condition_shape_shape(s);
            test     = LagrangianFunction.create(obj.boundaryMeshJoined, obj.meshDomain.ndim, 'P1'); % !!
            C        = lhs.compute(dLambda,test); 
        end

        function RHS = computeRHS(obj,u,dLambda)
%             F  = obj.computeForces(obj,k,u);
            F = zeros(u.nDofs,length(obj.DirFun));
            uD = obj.computeRHDcondition(dLambda,obj.DirFun);
            RHS = [F ; uD];
        end

        function RHS = computeForces(obj,stiffness,u)
            s.type      = 'Elastic';
            s.scale     = 'MACRO';
            s.dim.ndofs = u.nDofs;
            s.BC        = obj.boundaryConditions;
            s.mesh      = obj.meshDomain;
            RHSint      = RHSIntegrator.create(s);
            rhs         = RHSint.compute();
            % Perhaps move it inside RHSint?
            R           = RHSint.computeReactions(stiffness);
            RHS = rhs+R;
            RHS = obj.bcApplier.fullToReducedVectorDirichlet(RHS);
        end

        function uD = computeRHDcondition(obj,test,dir)
            s.mesh = obj.boundaryMeshJoined;
            s.quadType = 2;
            rhs = RHSIntegratorShapeFunction(s);
            nfun = size(dir,2);
            uD = [];
            for i=1:nfun
                rDiri = rhs.compute(dir{i},test);
                uD = [uD rDiri];
            end
        end

        function Meifem = createEIFEMPreconditioner(obj,mR,dir,iC,lG,bS,iCR,dMesh)
            % obj.EIFEMfilename = '/home/raul/Documents/Thesis/EIFEM/RAUL_rve_10_may_2024/EXAMPLE/EIFE_LIBRARY/DEF_Q4porL_2s_1.mat';
            EIFEMfilename = obj.fileNameEIFEM;
            % obj.EIFEMfilename = '/home/raul/Documents/Thesis/EIFEM/05_HEXAG2D/EIFE_LIBRARY/DEF_Q4auxL_1.mat';
            filename        = EIFEMfilename;
            s.RVE           = TrainedRVE(filename);
            s.mesh          = obj.createCoarseMesh(mR);
            s.DirCond       = dir;
            s.nSubdomains = obj.nSubdomains;
            eifem           = EIFEM(s);


            ss.ddDofManager = obj.createDomainDecompositionDofManager(iC,lG,bS,mR,iCR);
            ss.EIFEMsolver = eifem;
            ss.bcApplier = obj.bcApplier;
            ss.dMesh     = dMesh;
            ss.type = 'EIFEM';
            eP = Preconditioner.create(ss);
            Meifem = @(r,uk) eP.apply(r,uk);
        end

        function Mdn = createDirichletNeumannPreconditioner(obj,mR,dir,iC,lG,bS,lhs,mSb,iCR)
            s.ddDofManager  = obj.createDomainDecompositionDofManager(iC,lG,bS,mR,iCR);
            s.DirCond       = dir;
            s.bcApplier     = obj.bcApplier;
            s.LHS           = lhs;
            s.subdomainMesh = mSb;
            s.meshDomain = obj.meshDomain;
            s.type = 'DirichletNeumann';
            M = Preconditioner.create(s);
            Mdn = @(r) M.apply(r);
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

        function MgaussSeidel = createGaussSeidelpreconditioner(obj,LHS)
            s.LHS = LHS;
            s.type = 'GaussSeidel';
            M = Preconditioner.create(s);
            MgaussSeidel = @(r) M.apply(r);
        end

        function Mjacobi = createJacobipreconditioner(obj,LHS)
            s.LHS = LHS;
            s.type = 'Jacobi';
            M = Preconditioner.create(s);
            Mjacobi = @(r) M.apply(r);
        end


        function uD = AnalyticalDirCond(obj)
            test   = LagrangianFunction.create(obj.boundaryMeshJoined, obj.meshDomain.ndim, 'P1');
            s.mesh = obj.boundaryMeshJoined;
            s.quadType = 2;
%             rhs = RHSintegrator_ShapeFunction(s);

            f1x = @(x) [1/(4)*(1-x(1,:,:)).*(1-x(2,:,:));...
                    0*x(2,:,:)  ];
            f2x = @(x) [1/(4)*(1+x(1,:,:)).*(1-x(2,:,:));...
                    0*x(2,:,:)  ];
            f3x = @(x) [1/(4)*(1+x(1,:,:)).*(1+x(2,:,:));...
                    0*x(2,:,:)  ];
            f4x = @(x) [1/(4)*(1-x(1,:,:)).*(1+x(2,:,:));...
                    0*x(2,:,:)  ];

            f1y = @(x) [0*x(1,:,:);...
                    1/(4)*(1-x(1,:,:)).*(1-x(2,:,:))  ];
            f2y = @(x) [0*x(1,:,:);...
                    1/(4)*(1+x(1,:,:)).*(1-x(2,:,:))  ];
            f3y = @(x) [0*x(1,:,:);...
                    1/(4)*(1+x(1,:,:)).*(1+x(2,:,:))  ];
            f4y = @(x) [0*x(1,:,:);...
                    1/(4)*(1-x(1,:,:)).*(1+x(2,:,:))  ];

            f     = {f1x f1y f2x f2y f3x f3y f4x f4y}; %
            nfun = size(f,2);
            for i=1:nfun
                uD{i}  = AnalyticalFunction.create(f{i},obj.meshDomain.ndim,obj.boundaryMeshJoined);                   
            end


        end



%         function plotSolution(obj,x,mesh,row,col,iter,flag)
%             if nargin <7
%                 flag =0;
%             end
%             %             xFull = bc.reducedToFullVector(x);
%             if size(x,2)==1
%                 s.fValues = reshape(x,2,[])';
%             else
%                 s.fValues = x;
%             end
%             %
% 
%             s.mesh = mesh;
%             s.fValues(:,end+1) = 0;
%             s.ndimf = 2;
%             s.order = 'P1';
%             xF = LagrangianFunction(s);
%             %             xF.plot();
%             if flag == 0
%                 xF.print(['domain',num2str(row),num2str(col),'_',num2str(iter)],'Paraview')
%             elseif flag == 1
%                 xF.print(['DomainResidual',num2str(row),num2str(col),'_',num2str(iter)],'Paraview')
%             elseif flag == 2
%                 xF.print(['Residual',num2str(row),num2str(col),'_',num2str(iter)],'Paraview')
%             elseif flag == 3
%                 xF.print(['domainFine',num2str(row),num2str(col),'_',num2str(iter)],'Paraview')
%             elseif flag == 4
%                 xF.print(['domainNeuman',num2str(row),num2str(col),'_',num2str(iter)],'Paraview')
%             end
%             fclose('all');
%         end



    end

    methods (Static, Access = public)


        function J = computeTotalEnergy(x,A,b)
            J = 0.5*x'*A(x)-b'*x;
        end

        function plotSolution(x,mesh,row,col,iter,bcApplier,flag)
            if ~isempty(bcApplier)
                x = bcApplier.reducedToFullVectorDirichlet(x);
            end
            
            if nargin <7
                flag =0;
            end
            %             xFull = bc.reducedToFullVector(x);
            if size(x,2)==1
                s.fValues = reshape(x,2,[])';
            else
                s.fValues = x;
            end
            %

            s.mesh = mesh;
            s.fValues(:,end+1) = 0;
            s.ndimf = 2;
            s.order = 'P1';
            xF = LagrangianFunction(s);
            %             xF.plot();
            if flag == 0
                xF.print(['domain',num2str(row),num2str(col),'_',num2str(iter)],'Paraview')
            elseif flag == 1
                xF.print(['DomainResidual',num2str(row),num2str(col),'_',num2str(iter)],'Paraview')
            elseif flag == 2
                xF.print(['Residual',num2str(row),num2str(col),'_',num2str(iter)],'Paraview')
            elseif flag == 3
                xF.print(['domainFine',num2str(row),num2str(col),'_',num2str(iter)],'Paraview')
            elseif flag == 4
                xF.print(['domainNeuman',num2str(row),num2str(col),'_',num2str(iter)],'Paraview')
            end
            fclose('all');
        end
    end

end
