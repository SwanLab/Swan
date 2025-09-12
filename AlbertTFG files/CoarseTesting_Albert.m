classdef CoarseTesting_Albert < handle

    properties (Access = public)

    end
    properties (Access = private)
        nSubdomains
        rSubdomains
        mSubdomains
        ic
        icr
        lg
        bs

        meshDomain
        cellMeshes
        boundaryConditions
        bcApplier
        LHS
        RHS
        r
        centroids

        fileNameCorase
        tolSameNode

    end


    methods (Access = public)

        function obj = CoarseTesting_Albert()
            close all
            obj.init()

            mR = obj.createReferenceMesh();
            bS  = mR.createBoundaryMesh();
            [mD,mSb,iC,lG,iCR,discMesh] = obj.createMeshDomain(mR);
            obj.meshDomain = mD;
            obj.cellMeshes = mSb;
            obj.ic = iC;
            obj.icr = iCR;
            obj.lg = lG;
            obj.bs;

            [bC,dir] = obj.createBoundaryConditions(obj.meshDomain);
            obj.boundaryConditions = bC;
            obj.createBCapplier()

            [LHS,RHS,LHSf] = obj.createElasticProblem();
            obj.LHS = LHSf;
            %             LHS = 0.5*(LHS+LHS');

            LHSf = @(x) LHS*x;
            RHSf = RHS;
            Usol = LHS\RHS;
            Ufull = obj.bcApplier.reducedToFullVectorDirichlet(Usol);
            %obj.plotSolution(Ufull,obj.meshDomain,1,1,0,obj.bcApplier,0)

            % Meifem       = obj.createEIFEMPreconditioner(mR,dir,iC,lG,bS,iCR,discMesh);
            Milu         = obj.createILUpreconditioner(LHS);
            Mcoarse       = obj.createCoarsePreconditioner(mR,dir,iC,lG,bS,iCR,discMesh);
            Mid            = @(r) r;

            MiluCG = @(r,iter) Preconditioner.InexactCG(r,LHSf,Milu,RHSf);

            tol = 1e-8;
            tic
            x0 = zeros(size(RHSf));

            Mmult = @(r) Preconditioner.multiplePrec(r,Milu,Mcoarse,Milu,LHSf,RHSf,obj.meshDomain,obj.bcApplier);
            tic
            %           tau = @(r,A) 1;
            [uPCG,residualPCG,errPCG,errAnormPCG] = PCG.solve(LHSf,RHSf,x0,Mmult,tol,Usol,obj.meshDomain,obj.bcApplier);
            %            [uCG,residualPCG,errPCG,errAnormPCG] = RichardsonSolver.solve(LHSf,RHSf,x0,Mmult,tol,tau,Usol);
            toc
            
            xFull = obj.bcApplier.reducedToFullVectorDirichlet(uPCG);
            s.mesh = obj.meshDomain;
            s.ndimf = obj.meshDomain.ndim
            s.order = 'P1';
            s.fValues = reshape(xFull,2,[])';
            uFun = LagrangianFunction(s);

            obj.computeSubdomainCentroid();
            CoarsePlotSolution(uFun, obj.meshDomain, obj.bcApplier,'085 085 085', obj.r, obj.centroids);

            close all % Ho he afegit, cal borrar
            figure
            plot(residualPCG,'linewidth',2)
         %   hold on
         %   plot(residualCG,'linewidth',2)  SORTIA ERROR
         %   legend({'CG + ILU-EIFEM-ILU','CG'},'FontSize',12)
            set(gca, 'YScale', 'log')
            legend({'CG + ILU-EIFEM-ILU'},'FontSize',12)
            xlabel('Iteration')
            ylabel('Residual')

            figure
            plot(errPCG,'linewidth',2)
          %  hold on
          %  plot(errCG,'linewidth',2) SORTIA ERROR
            set(gca, 'YScale', 'log')
          %  legend('CG + EIFEM+ ILU(CG-90%-L2)','CG')
            legend('CG + EIFEM+ ILU(CG-90%-L2)')
            xlabel('Iteration')
            ylabel('||error||_{L2}')

            figure
            plot(errAnormPCG,'linewidth',2)
         %   hold on
         %   plot(errAnormCG,'linewidth',2)
            set(gca, 'YScale', 'log')
         %   legend('CG + EIFEM+ ILU(CG-90%-L2)','CG')
            legend('CG + EIFEM+ ILU(CG-90%-L2)')
            xlabel('Iteration')
            ylabel('Energy norm')

        end

    end

    methods (Access = private)

        function init(obj)
           % obj.nSubdomains    = [2 1]; %nx ny
            %obj.fileNameCorase = ["UL_r0_1-20x20.mat", "UL_r0_15-20x20.mat", "UL_r0_2-20x20.mat", "UL_r0_25-20x20.mat", "UL_r0_3-20x20.mat", "UL_r0_35-20x20.mat", "UL_r0_45-20x20.mat", "UL_r0_5-20x20.mat", "UL_r0_55-20x20.mat", "UL_r0_6-20x20.mat", "UL_r0_65-20x20.mat", "UL_r0_7-20x20.mat", "UL_r0_75-20x20.mat", "UL_r0_8-20x20.mat", "UL_r0_85-20x20.mat"]; %
            %obj.fileNameCorase = ["UL_r0_1-20x20.mat", "UL_r0_1-20x20.mat", "UL_r0_1-20x20.mat"];
            obj.fileNameCorase = ["UL_r0_85-20x20.mat", "UL_r0_85-20x20.mat", "UL_r0_85-20x20.mat"];
            %obj.fileNameCorase = ["UL_r0_85-20x20.mat","UL_r0_85-20x20.mat","UL_r0_85-20x20.mat","UL_r0_85-20x20.mat","UL_r0_85-20x20.mat","UL_r0_85-20x20.mat","UL_r0_85-20x20.mat","UL_r0_85-20x20.mat","UL_r0_85-20x20.mat","UL_r0_85-20x20.mat","UL_r0_85-20x20.mat","UL_r0_85-20x20.mat","UL_r0_85-20x20.mat","UL_r0_85-20x20.mat","UL_r0_85-20x20.mat"];
            %obj.fileNameCorase = ["UL_r0_1-20x20.mat","UL_r0_1-20x20.mat","UL_r0_1-20x20.mat","UL_r0_1-20x20.mat","UL_r0_1-20x20.mat","UL_r0_1-20x20.mat","UL_r0_1-20x20.mat","UL_r0_1-20x20.mat","UL_r0_1-20x20.mat","UL_r0_1-20x20.mat","UL_r0_1-20x20.mat","UL_r0_1-20x20.mat","UL_r0_1-20x20.mat","UL_r0_1-20x20.mat","UL_r0_85-20x20.mat"];
            %obj.fileNameCorase = ["UL_r0_1-20x20.mat","UL_r0_1-20x20.mat","UL_r0_85-20x20.mat"]
            %obj.rSubdomains    = [0.1, 0.1, 0.1, 0.1];
            obj.nSubdomains    = size(obj.fileNameCorase');
            obj.mSubdomains    = [];
            obj.tolSameNode    = 1e-10;
            obj.r = obj.loadRadius();
                      


        end

        function r = loadRadius(obj)
            r = zeros(size(obj.nSubdomains));
            for i = 1:obj.nSubdomains(1,2)
                for j = 1:obj.nSubdomains(1,1)
                    Data = load(obj.fileNameCorase(i,j));
                    r(i,j) = Data.R;
                end

            end

        end


        function [mD,mSb,iC,lG,iCR,discMesh] = createMeshDomain(obj,mR)
            s.nsubdomains   = obj.nSubdomains; %nx ny
            s.meshReference = mR;
            s.tolSameNode = obj.tolSameNode;
            m = MeshCreatorFromRVE(s);
            [mD,mSb,iC,~,lG,iCR,discMesh] = m.create();
        end


        function mS = createReferenceMesh(obj)
            mS = obj.createStructuredMesh();
            %   mS = obj.createMeshFromGid();
            % mS = obj.createEIFEMreferenceMesh();
        end


        function mS = createMeshFromGid(obj)
            filename   = 'lattice_ex1';
            a.fileName = filename;
            femD       = FemDataContainer(a);
            mS         = femD.mesh;
        end

        function mS = createStructuredMesh(obj)
             %UnitMesh better
            x1      = linspace(-1,1,20);
            x2      = linspace(-1,1,20);
            [xv,yv] = meshgrid(x1,x2);
            [F,V]   = mesh2tri(xv,yv,zeros(size(xv)),'x');
            s.coord  = V(:,1:2);
            s.connec = F;
            mS = Mesh.create(s);
          
        end

        function computeSubdomainCentroid(obj)
            for i = 1:obj.nSubdomains(1,2)
                for j = 1:obj.nSubdomains(1,1)
                   x0=mean(obj.cellMeshes{i,j}.coord(:,1));
                   y0=mean(obj.cellMeshes{i,j}.coord(:,2));
                   obj.centroids = cat(1,obj.centroids, [x0,y0]);
                end

            end

            

        end

        function levelSet = createLevelSetFunction(obj,bgMesh)
            sLS.type        = 'CircleInclusion';
            sLS.xCoorCenter = 0;
            sLS.yCoorCenter = 0;
            sLS.radius      = 0.1;
            g               = GeometricalFunction(sLS);
            lsFun           = g.computeLevelSetFunction(bgMesh);
            levelSet        = lsFun.fValues;
        end
        
        function uMesh = computeUnfittedMesh(~,bgMesh,levelSet)
            sUm.backgroundMesh = bgMesh;
            sUm.boundaryMesh   = bgMesh.createBoundaryMesh();
            uMesh              = UnfittedMesh(sUm);
            uMesh.compute(levelSet);
        end


        function mS = createEIFEMreferenceMesh(obj)
            filename = obj.fileNameCorase;
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
            % connec = [2 3 4 1];
            s.coord = coord;
            s.connec = connec;
            cMesh = Mesh.create(s);
        end

        function createBCapplier(obj)
            s.mesh                  = obj.meshDomain;
            s.boundaryConditions    = obj.boundaryConditions;
            obj.bcApplier           = BCApplier(s);
        end


        function material = createMaterial(obj)
            
            material = cell(size(obj.fileNameCorase));

            for i = 1:obj.nSubdomains(1,2)
                for j = 1:obj.nSubdomains(1,1)
                    [young,poisson] = obj.computeElasticProperties(obj.cellMeshes{i,j}, obj.r(i,j) );

                    s.type        = 'ISOTROPIC';
                    s.ptype       = 'ELASTIC';
                    s.ndim        = obj.cellMeshes{i,j}.ndim;
                    s.young       = young;
                    s.poisson     = poisson;
                    tensor        = Material.create(s);
                    material{i,j} = tensor;

                end
            end
      

        end


        function material = createMaterialBasic(obj,mesh)
            [young,poisson] = obj.computeElasticProperties(mesh);
            s.type    = 'ISOTROPIC';
            s.ptype   = 'ELASTIC';
            s.ndim    = mesh.ndim;
            s.young   = young;
            s.poisson = poisson;
            tensor    = Material.create(s);
            material  = tensor;
        end

        function [young,poisson] = computeElasticProperties(obj,mesh, radius)
            E1  = 1;
            E2 = E1/1000;
            nu = 1/3; 
            x0=mean(mesh.coord(:,1));
            y0=mean(mesh.coord(:,2));
%             young   = ConstantFunction.create(E,mesh);
%             poisson = ConstantFunction.create(nu,mesh);
            f   = @(x) (sqrt((x(1,:,:)-x0).^2+(x(2,:,:)-y0).^2)<radius)*E2 + ...
                        (sqrt((x(1,:,:)-x0).^2+(x(2,:,:)-y0).^2)>=radius)*E1 ; 
%                                      x(2,:,:).*0 ];
            young   = AnalyticalFunction.create(f,1,mesh);
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
                        PL.direction = [2];
                        PL.value     = [1];       %Set displacement intensity ------------------------------------------------------------
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
            u = LagrangianFunction.create(obj.cellMeshes{1,1},obj.cellMeshes{1,1}.ndim,'P1');
            uBasic = LagrangianFunction.create(obj.meshDomain,obj.meshDomain.ndim,'P1');
            material = obj.createMaterial();
            %materialBasic = obj.createMaterialBasic(obj.meshDomain);
            [lhs,LHSr] = obj.computeStiffnessMatrix(u,material);
            %[lhsBasic,~] = obj.computeStiffnessMatrixBasic(uBasic,materialBasic);
            RHSr       = obj.computeForces(lhs,uBasic);
        end


        function [LHS,LHSr] = computeStiffnessMatrix(obj,dispFun,mat)
            s.type     = 'ElasticStiffnessMatrix';
            s.quadratureOrder = 2;
            s.test     = dispFun;
            s.trial    = dispFun;
            % LHScell    = cell(size(obj.rSubdomains));
            % LHScellGlobal = LHScell;
            LHSvect = [];

            for i = 1:obj.nSubdomains(1,2)
                for j = 1:obj.nSubdomains(1,1)
                    s.mesh     = obj.cellMeshes{i,j};
                    s.material = mat{i,j};
                    
                    lhs          = LHSIntegrator.create(s);
                    % LHScell{i,j} = full(lhs.compute());
                    LHSvect = cat(3, LHSvect, full(lhs.compute()) );
                    
                end
            end

            p.nSubdomains = obj.nSubdomains;
            p.interfaceConnec = obj.ic;
            p.interfaceConnecReshaped = obj.icr;
            p.locGlobConnec = obj.lg;
            p.nBoundaryNodes = obj.bs;
            p.nReferenceNodes = obj.cellMeshes{1,1}.nnodes;
            p.nNodes = obj.meshDomain.nnodes;
            p.nDimf = obj.meshDomain.ndim;
            
            
            
            dddm = DomainDecompositionDofManager(p);
            LHS = dddm.local2globalMatrix(LHSvect);
            LHSr = obj.bcApplier.fullToReducedMatrixDirichlet(LHS);
            
        end


        function [LHS,LHSr] = computeStiffnessMatrixBasic(obj,dispFun,mat)
            s.type     = 'ElasticStiffnessMatrix';
            s.mesh     = obj.meshDomain;
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
            EIFEMfilename = obj.fileNameCorase;
            % obj.EIFEMfilename = '/home/raul/Documents/Thesis/EIFEM/05_HEXAG2D/EIFE_LIBRARY/DEF_Q4auxL_1.mat';
            filename        = EIFEMfilename;
            s.RVE           = CoarseTrainedRVE(filename);
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

         function Mcoarse = createCoarsePreconditioner(obj,mR,dir,iC,lG,bS,iCR,dMesh)
            % obj.EIFEMfilename = '/home/raul/Documents/Thesis/EIFEM/RAUL_rve_10_may_2024/EXAMPLE/EIFE_LIBRARY/DEF_Q4porL_2s_1.mat';
            % obj.EIFEMfilename = '/home/raul/Documents/Thesis/EIFEM/05_HEXAG2D/EIFE_LIBRARY/DEF_Q4auxL_1.mat';
            
            RVE = cell(size(obj.nSubdomains));

            for i = 1:obj.nSubdomains(1,2)
                for j = 1:obj.nSubdomains(1,1)
                    filename = obj.fileNameCorase(i,j);
                    RVE{i,j} = CoarseTrainedRVE(filename);  %%Passar vector de filenames
                end
            end
            
            
            s.RVE           = RVE;
            s.mesh          = obj.createCoarseMesh(obj.cellMeshes{1,1});
            s.DirCond       = dir;
            s.nSubdomains   = obj.nSubdomains;
            coarseSolver    = Coarse(s);


            ss.ddDofManager = obj.createDomainDecompositionDofManager(iC,lG,bS,mR,iCR);
            ss.Coarsesolver = coarseSolver;
            ss.bcApplier = obj.bcApplier;
            ss.dMesh     = dMesh;
            ss.type = 'Coarse';
            eP = Preconditioner.create(ss);
            Mcoarse = @(r) eP.apply(r);
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
            s.LHS = sparse(LHS);
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
