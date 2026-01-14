classdef EIFEMtesting_3D < handle

    properties (Access = public)

    end
    properties (Access = private)
        meshDomain
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

        function obj = EIFEMtesting_3D()
            close all
            obj.init()

            mR = obj.createReferenceMesh();
            bS  = mR.createBoundaryMesh();
            [mD,mSb,iC,lG,iCR,discMesh] = obj.createMeshDomain(mR);
            obj.meshDomain = mD;
            [bC,dir] = obj.createBoundaryConditions(obj.meshDomain);
            obj.boundaryConditions = bC;
            obj.createBCapplier()

            tic
            [LHS,RHS,LHSf] = obj.createElasticProblem();
            toc
            obj.LHS = LHSf;
            %             LHS = 0.5*(LHS+LHS');

            LHSf = @(x) LHS*x;
            RHSf = RHS;
            rhs2 = repmat(RHS,[1,24,1]);
            tic
            Usol = LHS\RHS;
            toc
            Ufull = obj.bcApplier.reducedToFullVectorDirichlet(Usol);

            %obj.plotSolution(Ufull,obj.meshDomain,1,1,0,obj.bcApplier,0)

            RBbasisFree  = forAlgebraicMultigrid(obj);
            Mid          = @(r) r;
            Meifem       = obj.createEIFEMPreconditioner(mR,dir,iC,lG,bS,iCR,discMesh);
            %Milu         = obj.createILUpreconditioner(LHS);
            MgaussSeidel = obj.createGaussSeidelpreconditioner(LHS);
            MJacobi      = obj.createJacobipreconditioner(LHS);
            Mmodal       = obj.createModalpreconditioner(LHS);
            %            MblockD      = obj.createBlockDiagonalpreconditioner(LHS);
            %             MdirNeu      = obj.createDirichletNeumannPreconditioner(mR,dir,iC,lG,bS,obj.LHS,mSb);

            MiluCG = @(r,iter) Preconditioner.InexactCG(r,LHSf,Milu,RHSf);

            tol = 1e-8;

            x0 = zeros(size(RHSf));
            tic
            [uCG,residualCG,errCG,errAnormCG] = PCG.solve(LHSf,RHSf,x0,MJacobi ...
                ,tol,Usol,obj.meshDomain,obj.bcApplier);
            toc
            %             [uCG,residualCG,errCG,errAnormCG] = RichardsonSolver.solve(LHSf,RHSf,x0,P,tol,0.1,Usol);

            tol = 1e-6;

            %Mmult = MdirNeu;
            x0 = zeros(size(RHSf));
            r = RHSf - LHSf(x0);
            Mmult = @(r) Preconditioner.multiplePrec(r,Mid,Meifem,Mid,LHSf,RHSf,obj.meshDomain,obj.bcApplier);
            %              Mmult = @(r) Preconditioner.multiplePrec(r,Mid,Meifem,Mid,LHSf,RHSf,obj.meshDomain,obj.bcApplier);
            %             zmult = Mmult(r);

            %             zfull = obj.bcApplier.reducedToFullVectorDirichlet(zmult);
            %obj.plotSolution(zfull,obj.meshDomaopenin,0,0,2,obj.bcApplier,0)

            %             zeifem = Meifem(r);
            %             zfull = obj.bcApplier.reducedToFullVectorDirichlet(zeifem);
            %obj.plotSolution(zfull,obj.meshDomain,0,0,1,obj.bcApplier,0)
            % x0 = zmult;
            tic
            %           tau = @(r,A) 1;
            [uPCG,residualPCG,errPCG,errAnormPCG] = PCG.solve(LHSf,RHSf,x0,Mmult,tol,Usol,obj.meshDomain,obj.bcApplier);
            %            [uCG,residualPCG,errPCG,errAnormPCG] = RichardsonSolver.solve(LHSf,RHSf,x0,Mmult,tol,tau,Usol);
            toc

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

    methods (Access = private)

        function init(obj)
            obj.nSubdomains  = [5 1 1]; %nx ny
            obj.fileNameEIFEM = 'MESHdom.mat';
            %             obj.fileNameEIFEM = 'DEF_auxNew_2.mat';
            %obj.fileNameEIFEM = 'DEF_Q4porL_1_raul.mat';
            obj.tolSameNode = 1e-6;
        end

        function [mD,mSb,iC,lG,iCR,discMesh] = createMeshDomain(obj,mR)
            s.nsubdomains   = obj.nSubdomains; %nx ny
            s.meshReference = mR;
            s.tolSameNode = obj.tolSameNode;
            m = MeshCreatorFromRVE3D(s);
            [mD,mSb,iC,~,lG,iCR,discMesh] = m.create();
        end


        function mS = createReferenceMesh(obj)
            %                             mS = obj.createStructuredMesh();
            %             mS = obj.createMeshFromGid();
            mS = obj.createEIFEMreferenceMesh();
        end


        function mS = createMeshFromGid(obj)
            %             filename   = 'lattice_ex1';
            %             a.fileName = filename;
            %             femD       = FemDataContainer(a);
            %             mS         = femD.mesh;
            load('por3Dmesh.mat');
            s.connec = porMesh.connec ;
            s.coord  = porMesh.coord;
            maxC= max(s.coord);
            minC = min(s.coord);
            s.coord(s.coord(:,1)== maxC(1) & s.coord(:,3)==maxC(3),:) =...
                s.coord(s.coord(:,1)== maxC(1) & s.coord(:,3)==maxC(3),:)-[0,0,1e-5];

            s.coord(s.coord(:,1)== maxC(1) & s.coord(:,3)==minC(3),:) =...
                s.coord(s.coord(:,1)== maxC(1) & s.coord(:,3)==minC(3),:)+[0,0,1e-5];

            s.coord(s.coord(:,1)== minC(1) & s.coord(:,3)==maxC(3),:) =...
                s.coord(s.coord(:,1)== minC(1) & s.coord(:,3)==maxC(3),:)-[0,0,1e-5];

            s.coord(s.coord(:,1)== minC(1) & s.coord(:,3)==minC(3),:) =...
                s.coord(s.coord(:,1)== minC(1) & s.coord(:,3)==minC(3),:)+[0,0,1e-5];

            s.coord(s.coord(:,2)== maxC(2) & s.coord(:,3)==maxC(3),:) =...
                s.coord(s.coord(:,2)== maxC(2) & s.coord(:,3)==maxC(3),:)-[0,0,1e-5];

            s.coord(s.coord(:,2)== maxC(2) & s.coord(:,3)==minC(3),:) =...
                s.coord(s.coord(:,2)== maxC(2) & s.coord(:,3)==minC(3),:)+[0,0,1e-5];

            s.coord(s.coord(:,2)== minC(2) & s.coord(:,3)==maxC(3),:) =...
                s.coord(s.coord(:,2)== minC(2) & s.coord(:,3)==maxC(3),:)-[0,0,1e-5];

            s.coord(s.coord(:,2)== minC(2) & s.coord(:,3)==minC(3),:) =...
                s.coord(s.coord(:,2)== minC(2) & s.coord(:,3)==minC(3),:)+[0,0,1e-5];

            mS = Mesh.create(s);


        end

        function mS = createStructuredMesh(obj)
            %             mS = HexaMesh(1,1,1,10,10,10);
            mS = TetraMesh(1,1,1,10,10,10);
            s.coord = mS.coord;
            maxC= max(s.coord);
            minC = min(s.coord);
            s.coord(s.coord(:,1)== maxC(1) & s.coord(:,3)==maxC(3),:) =...
                s.coord(s.coord(:,1)== maxC(1) & s.coord(:,3)==maxC(3),:)-[0,0,1e-5];

            s.coord(s.coord(:,1)== maxC(1) & s.coord(:,3)==minC(3),:) =...
                s.coord(s.coord(:,1)== maxC(1) & s.coord(:,3)==minC(3),:)+[0,0,1e-5];

            s.coord(s.coord(:,1)== minC(1) & s.coord(:,3)==maxC(3),:) =...
                s.coord(s.coord(:,1)== minC(1) & s.coord(:,3)==maxC(3),:)-[0,0,1e-5];

            s.coord(s.coord(:,1)== minC(1) & s.coord(:,3)==minC(3),:) =...
                s.coord(s.coord(:,1)== minC(1) & s.coord(:,3)==minC(3),:)+[0,0,1e-5];

            s.coord(s.coord(:,2)== maxC(2) & s.coord(:,3)==maxC(3),:) =...
                s.coord(s.coord(:,2)== maxC(2) & s.coord(:,3)==maxC(3),:)-[0,0,1e-5];

            s.coord(s.coord(:,2)== maxC(2) & s.coord(:,3)==minC(3),:) =...
                s.coord(s.coord(:,2)== maxC(2) & s.coord(:,3)==minC(3),:)+[0,0,1e-5];

            s.coord(s.coord(:,2)== minC(2) & s.coord(:,3)==maxC(3),:) =...
                s.coord(s.coord(:,2)== minC(2) & s.coord(:,3)==maxC(3),:)-[0,0,1e-5];

            s.coord(s.coord(:,2)== minC(2) & s.coord(:,3)==minC(3),:) =...
                s.coord(s.coord(:,2)== minC(2) & s.coord(:,3)==minC(3),:)+[0,0,1e-5];
            s.connec = mS.connec;
            mS = Mesh.create(s);
            u = LagrangianFunction.create(mS,mS.ndim,'P1');
            %             % Generate coordinates
            %             x1 = linspace(0,1,2);
            %             x2 = linspace(0,1,2);
            %             % Create the grid
            %             [xv,yv] = meshgrid(x1,x2);
            %             % Triangulate the mesh to obtain coordinates and connectivities
            %             [F,coord] = mesh2tri(xv,yv,zeros(size(xv)),'x');
            %
            %             s.coord    = coord(:,1:2);
            %             s.connec   = F;
            %             mS         = Mesh.create(s);
        end

        function mS = createEIFEMreferenceMesh(obj)
            filename = obj.fileNameEIFEM;
            load(filename);
            s.coord    = MESHdom.COOR;
            s.connec   = MESHdom.CN;

            maxC= max(s.coord);
            minC = min(s.coord);
            s.coord(s.coord(:,1)== maxC(1) & s.coord(:,3)==maxC(3),:) =...
                s.coord(s.coord(:,1)== maxC(1) & s.coord(:,3)==maxC(3),:)-[0,0,1e-5];

            s.coord(s.coord(:,1)== maxC(1) & s.coord(:,3)==minC(3),:) =...
                s.coord(s.coord(:,1)== maxC(1) & s.coord(:,3)==minC(3),:)+[0,0,1e-5];

            s.coord(s.coord(:,1)== minC(1) & s.coord(:,3)==maxC(3),:) =...
                s.coord(s.coord(:,1)== minC(1) & s.coord(:,3)==maxC(3),:)-[0,0,1e-5];

            s.coord(s.coord(:,1)== minC(1) & s.coord(:,3)==minC(3),:) =...
                s.coord(s.coord(:,1)== minC(1) & s.coord(:,3)==minC(3),:)+[0,0,1e-5];

            s.coord(s.coord(:,2)== maxC(2) & s.coord(:,3)==maxC(3),:) =...
                s.coord(s.coord(:,2)== maxC(2) & s.coord(:,3)==maxC(3),:)-[0,0,1e-5];

            s.coord(s.coord(:,2)== maxC(2) & s.coord(:,3)==minC(3),:) =...
                s.coord(s.coord(:,2)== maxC(2) & s.coord(:,3)==minC(3),:)+[0,0,1e-5];

            s.coord(s.coord(:,2)== minC(2) & s.coord(:,3)==maxC(3),:) =...
                s.coord(s.coord(:,2)== minC(2) & s.coord(:,3)==maxC(3),:)-[0,0,1e-5];

            s.coord(s.coord(:,2)== minC(2) & s.coord(:,3)==minC(3),:) =...
                s.coord(s.coord(:,2)== minC(2) & s.coord(:,3)==minC(3),:)+[0,0,1e-5];

            mS = Mesh.create(s);


        end

        function mCoarse = createCoarseMesh(obj,mR)
            s.nsubdomains   = obj.nSubdomains; %nx ny
            s.meshReference = obj.createReferenceCoarseMesh(mR);
            %             s.meshReference = obj.loadReferenceCoarseMesh(mR);
            s.tolSameNode   = obj.tolSameNode;
            mRVECoarse      = MeshCreatorFromRVE3D(s);
            [mCoarse,~,~] = mRVECoarse.create();
        end



        function cMesh = createReferenceCoarseMesh(obj,mR)
            xmax = max(mR.coord(:,1));
            xmin = min(mR.coord(:,1));
            ymax = max(mR.coord(:,2));
            ymin = min(mR.coord(:,2));
            zmax = max(mR.coord(:,3));
            zmin = min(mR.coord(:,3));
            coord(1,1) = xmax;  coord(1,2) = ymax;   coord(1,3) = zmax;
            coord(2,1) = xmin;  coord(2,2) = ymax;   coord(2,3) = zmax;
            coord(3,1) = xmin;  coord(3,2) = ymax;   coord(3,3) = zmin;
            coord(4,1) = xmax;  coord(4,2) = ymax;   coord(4,3) = zmin;
            coord(5,1) = xmax;  coord(5,2) = ymin;   coord(5,3) = zmax;
            coord(6,1) = xmin;  coord(6,2) = ymin;   coord(6,3) = zmax;
            coord(7,1) = xmin;  coord(7,2) = ymin;   coord(7,3) = zmin;
            coord(8,1) = xmax;  coord(8,2) = ymin;   coord(8,3) = zmin;

            connec = [1 2 3 4 5 6 7 8];
            s.coord = coord;
            s.connec = connec;
            cMesh = Mesh.create(s);
        end

        function cMesh = loadReferenceCoarseMesh(obj,mR)
            bS  = mR.createBoundaryMesh();
            bS2{1} = bS{3}; bS2{2} = bS{2}; bS2{3} = bS{4}; bS2{4} = bS{1}; % reorder boundaries
            bS = bS2;
            nbd = size(bS,2);
            interpType = [2,1,2,1];
            inode = 1;
            for ibd = 1:nbd
                maxCoord  = max(bS{ibd}.mesh.coord);
                minCoord  = min(bS{ibd}.mesh.coord);
                meanCoord = (maxCoord+minCoord)/2;
                val = ibd<=nbd/2;
                if interpType(ibd) == 1
                    coord(inode,:)   = val*minCoord + abs((val-1))*maxCoord;
                    coord(inode+1,:) = val*maxCoord + abs((val-1))*minCoord;
                    inode=inode + 2;
                else
                    coord(inode,:)   = val*minCoord + abs((val-1))*maxCoord;
                    coord(inode+1,:) = meanCoord;
                    coord(inode+2,:) = val*maxCoord + abs((val-1))*minCoord;
                    inode=inode + 3;
                end
            end

            %              coord(1,:)  = [ 0.378041543026706 , -0.843442136498517 ];
            %              coord(2,:)  = [ 1.49050445103858  , -0.843442136498517 ];
            %              coord(3,:)  = [ 2.60296735905045  , -0.843442136498517 ];
            %              coord(4,:)  = [ 2.98100890207715  ,  0                 ];
            %              coord(5,:)  = [ 2.98100890207715  ,  0.314540059347181 ];
            %              coord(6,:)  = [ 2.60296735905045  ,  1.1579821958457   ];
            %              coord(7,:)  = [ 1.49050445103858  ,  1.1579821958457   ];
            %              coord(8,:)  = [ 0.378041543026706 ,  1.1579821958457   ];
            %              coord(9,:)  = [ 0                 ,  0.314540059347181 ];
            %              coord(10,:) = [ 0                 ,  0                 ];
            %
            %


            connec = [1 2 3 4 5 6 7 8 9 10];
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
            %             E  = 70000;
            %             nu = 0.3;
            Epstr  = E/(1-nu^2);
            nupstr = nu/(1-nu);
            %             young   = ConstantFunction.create(Epstr,mesh);
            %             poisson = ConstantFunction.create(nupstr,mesh);
            young   = ConstantFunction.create(E,mesh);
            poisson = ConstantFunction.create(nu,mesh);
        end

        function [Dir,PL] = createRawBoundaryConditions(obj)
            minx = min(obj.meshDomain.coord(:,1));
            maxx = max(obj.meshDomain.coord(:,1));
            miny = min(obj.meshDomain.coord(:,2));
            maxy = max(obj.meshDomain.coord(:,2));
            minz = min(obj.meshDomain.coord(:,3));
            maxz = max(obj.meshDomain.coord(:,3));
            tolBound = obj.tolSameNode;
            isLeft   = @(coor) (abs(coor(:,1) - minx)   < tolBound);
            isRight  = @(coor) (abs(coor(:,1) - maxx)   < tolBound);
            isBottom = @(coor) (abs(coor(:,3) - minz)   < tolBound);
            isTop    = @(coor) (abs(coor(:,3) - maxz)   < tolBound);
            %             isMiddle = @(coor) (abs(coor(:,2) - max(coor(:,2)/2)) == 0);
            Dir{1}.domain    = @(coor) isLeft(coor);%| isRight(coor) ;
            Dir{1}.direction = [1,2,3];
            Dir{1}.value     = 0;

            %             Dir{2}.domain    = @(coor) isRight(coor) ;
            %             Dir{2}.direction = [2];
            %             Dir{2}.value     = 0;

            %             PL.domain    = @(coor) isTop(coor);
            %             PL.direction = [2];
            %             PL.value     = [-0.1];
            PL.domain    = @(coor) isRight(coor);
            PL.direction = [2];
            PL.value     = [-1];
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


        function [LHSr,RHSr,lhs] = createElasticProblem(obj)
            u = LagrangianFunction.create(obj.meshDomain,obj.meshDomain.ndim,'P1');
            material = obj.createMaterial(obj.meshDomain);
            [lhs,LHSr] = obj.computeStiffnessMatrix(obj.meshDomain,u,material);
            RHSr       = obj.computeForces(lhs,u);
        end


        function [LHS,LHSr] = computeStiffnessMatrix(obj,mesh,dispFun,mat)
            s.type     = 'ElasticStiffnessMatrix';
            s.mesh     = mesh;
            s.test     = dispFun;
            s.trial    = dispFun;
            s.material = mat;
            s.quadratureOrder = 2;
            lhs = LHSIntegrator.create(s);
            LHS = lhs.compute();
            LHSr = obj.bcApplier.fullToReducedMatrixDirichlet(LHS);
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

        function Meifem = createEIFEMPreconditioner(obj,mR,dir,iC,lG,bS,iCR,dMesh)
            % obj.EIFEMfilename = '/home/raul/Documents/Thesis/EIFEM/RAUL_rve_10_may_2024/EXAMPLE/EIFE_LIBRARY/DEF_Q4porL_2s_1.mat';
            EIFEMfilename = 'DEF_Q8_wing_1.mat';
            % obj.EIFEMfilename = '/home/raul/Documents/Thesis/EIFEM/05_HEXAG2D/EIFE_LIBRARY/DEF_Q4auxL_1.mat';
            filename        = EIFEMfilename;
            s.RVE           = TrainedRVE(filename);
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
            d = DomainDecompositionDofManager3D(s);
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

        function Mmodal = createModalpreconditioner(obj,LHS)
            s.LHS = LHS;
            s.nBasis = 8;
            s.type   = 'MODAL';
            M = Preconditioner.create(s);
            Mmodal = @(r) M.apply(r);
        end

        function MblockD = createBlockDiagonalpreconditioner(obj,LHS)
            s.LHS       = LHS;
            s.dimension = 20;
            s.type      = 'BlockDiagonal';
            M = Preconditioner.create(s);
            MblockD = @(r) M.apply(r);
        end


        function B = forAlgebraicMultigrid(obj)
%             refPoint = (min(obj.meshDomain.coord)+max(obj.meshDomain.coord))/2 ;
             refPoint = min(obj.meshDomain.coord);
             Bfull = obj.rigidBodyModes3D(obj.meshDomain.coord,refPoint);
             for i = 1:size(Bfull,2)
                 B(:,i) = obj.bcApplier.fullToReducedVectorDirichlet(Bfull(:,i));
             end


% 
%             RB = RigidBodyFunction.create(obj.meshDomain,refPoint);
%             xt  = RB.basisFunctions{1}.project('P1');
%             yt  = RB.basisFunctions{2}.project('P1');
%             rot  = RB.basisFunctions{3}.project('P1');
%             BG = [reshape(xt.fValues',[],1),reshape(yt.fValues',[],1),reshape(rot.fValues',[],1)];
%             B = [obj.bcApplier.fullToReducedVectorDirichlet(BG(:,1)),obj.bcApplier.fullToReducedVectorDirichlet(BG(:,2)),...
%                 obj.bcApplier.fullToReducedVectorDirichlet(BG(:,3))];
        end

        function R = rigidBodyModes3D(obj,coords, refPoint)
            % rigidBodyModes3D computes the 6 rigid body modes in 3D
            %
            % INPUTS:
            %   coords    : N x 3 matrix of nodal coordinates [x, y, z]
            %   refPoint  : 1 x 3 vector [x0, y0, z0] - reference point for rotation
            %
            % OUTPUT:
            %   R         : 3N x 6 matrix where each column is a rigid body mode
            %               (3 displacements per node)

            % Number of nodes
            N = size(coords, 1);

            % Initialize rigid body modes matrix
            R = zeros(3*N, 6);

            % Displacement indices
            ix = 1:3:3*N;
            iy = 2:3:3*N;
            iz = 3:3:3*N;

            % Displacements for translation modes
            R(ix, 1) = 1; % Translation in x
            R(iy, 2) = 1; % Translation in y
            R(iz, 3) = 1; % Translation in z

            % Compute position vectors relative to reference point
            relCoords = coords - refPoint;

            x = relCoords(:, 1);
            y = relCoords(:, 2);
            z = relCoords(:, 3);

            % Displacements for rotation about x-axis
            R(iy, 4) = -z;
            R(iz, 4) =  y;

            % Displacements for rotation about y-axis
            R(ix, 5) =  z;
            R(iz, 5) = -x;

            % Displacements for rotation about z-axis
            R(ix, 6) = -y;
            R(iy, 6) =  x;

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
                s.fValues = reshape(x,mesh.ndim,[])';
            else
                s.fValues = x;
            end

            if mesh.ndim == 2
                s.fValues(:,end+1) = 0;
            end
            %

            s.mesh = mesh;
            s.ndimf = mesh.ndim;
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
