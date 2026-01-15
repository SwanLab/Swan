classdef ModalWingEIFEM3D < handle

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

        function obj = ModalWingEIFEM3D()
            close all
            obj.init()

            mR = obj.createReferenceMesh();
            bS  = mR.createBoundaryMesh();
            [mD,mSb,iC,lG,iCR,discMesh] = obj.createMeshDomain(mR);
            obj.meshDomain = mD;
            [bC,dir] = obj.createBoundaryConditions(obj.meshDomain);
            obj.boundaryConditions = bC;
            obj.createBCapplier();           

            [LHS, Mr, RHS, LHSf] = obj.createElasticProblem();

            LHSf = @(x) LHS*x;

            RHSf = RHS;

            Usol = LHS\RHS;
            Ufull = obj.bcApplier.reducedToFullVectorDirichlet(Usol);

            Meifem       = obj.createEIFEMPreconditioner(mR,dir,iC,lG,bS,iCR,discMesh);

            Milu         = obj.createILUpreconditioner((LHS'+LHS)/2);
            
             Mmult = @(r) Preconditioner.multiplePrec(r,Milu,Meifem,Milu,LHSf,RHSf,obj.meshDomain,obj.bcApplier);
            Mid          = @(r) r;
            %% Solve Stastic Problem

            
            tol = 1e-8;
            x0 = zeros(size(RHSf));
            tic
            [uCG,residualCG,errCG,errAnormCG] = PCG.solve(LHSf,RHS,x0,Mid,tol,Usol,obj.meshDomain,obj.bcApplier);
            toc

            x0 = zeros(size(RHSf));

            tic
            [uPCG,residualPCG,errPCG,errAnormPCG] = PCG.solve(LHSf,RHSf,x0,Mmult,tol,Usol,obj.meshDomain,obj.bcApplier);
            toc

            figure
            plot(residualPCG,'linewidth',2)
            hold on
            plot(residualCG,'linewidth',2)
            set(gca, 'YScale', 'log')
            legend({'CG + ILU-EIFEM-ILU','CG'},'FontSize',12)
            xlabel('Iteration')
            ylabel('Residual')


            %% Solver modal Problem


            problem = 3; %1 = 4DOF spring, 2 =) 1000DOF Spring, 3 = testingEIFEM matrices

            solver = LOBPCG(problem,LHS,Mr,Mmult);
            % Solver Controls:
            solver.b      = 5;        % number of modes 
            solver.use_precond = true;
            solver.use_physical_x = false; % feed aproximation to eigenvaluees (simulating coarse from superelement) as starting point
            solver.tol    = 1e-8;    % residual tolerance
            solver.maxit  = 20000;       % max iterations
            solver.verbose = false;
            solver.precond_type = 'eifem';

            results = solver.run_demo();

            %% Jacobi with EIFEM initial guess 
            solver2 = LOBPCG(problem,LHS,Mr);
            solver2.b      = 5;        % number of modes 
            solver2.use_precond = true;
            solver2.use_physical_x = true; % feed aproximation to eigenvaluees (simulating coarse from superelement) as starting point
            solver2.tol    = 1e-8;    % residual tolerance
            solver2.maxit  = 20000;       % max iterations
            solver2.verbose = false;
            solver2.precond_type ='jacobi';
            results2 = solver2.run_demo();

            tol = 1e-8; 
            tol2 = 1e-8;
            for i = 1:5
                figure(i)
                
                % --- Strategy 1: Jacobi + EIFEM Init
                raw_res1 = results2.history.rnorm(:,i);
                % Filter: keep only values strictly greater than tolerance
                res_1 = raw_res1(raw_res1 > tol); 
                x_axis1 = 1:length(res_1);
            
                % --- Strategy 2: No Prec + EIFEM Init (results4) ---
                raw_res2 = results.history.rnorm(:,i);
                res_2 = raw_res2(raw_res2 > tol2);
                x_axis2 = 1:length(res_2);
            
                % Plotting with dynamic lengths
                semilogy(x_axis1, res_1, ...
                         x_axis2, res_2, ...
                         'LineWidth', 1)
                         
                xlabel('Iterations')
                ylabel('Residual')
                grid on
                title(['Residual Evolution For Different Strategies (Mode ' num2str(i) ')'])
                
                legend('Jacobi + EIFEM Init', ...
                       'EIFEM + EIFEM Init', ...
                       'Location', 'best') 
            end
            
            obj.printFineSolution(results.X,mD);
            

        end

    end

    methods (Access = private)

        function init(obj)
            obj.nSubdomains  = [5 1 1]; %nx ny
            obj.fileNameEIFEM = 'DEF_Q8_wing_1.mat';
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
            mS = obj.createEIFEMreferenceMesh();
        end

        function mS = createEIFEMreferenceMesh(obj)
            filename = obj.fileNameEIFEM;
            load(filename);
            s.coord    = EIFEoper.MESH.COOR;
            s.connec   = EIFEoper.MESH.CN;

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
            s.nsubdomains   = obj.nSubdomains; %nx ny nz
            s.meshReference = obj.createReferenceCoarseMesh(mR);
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

            coord(1,1) = xmin;  coord(1,2) = ymin;   coord(1,3) = zmax;
            coord(2,1) = xmin;  coord(2,2) = ymax;   coord(2,3) = zmax;
            coord(3,1) = xmax;  coord(3,2) = ymax;   coord(3,3) = zmax;
            coord(4,1) = xmax;  coord(4,2) = ymin;   coord(4,3) = zmax;
            coord(5,1) = xmin;  coord(5,2) = ymin;   coord(5,3) = zmin;
            coord(6,1) = xmin;  coord(6,2) = ymax;   coord(6,3) = zmin;
            coord(7,1) = xmax;  coord(7,2) = ymax;   coord(7,3) = zmin;
            coord(8,1) = xmax;  coord(8,2) = ymin;   coord(8,3) = zmin;

            connec = [1 2 3 4 5 6 7 8];
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

            E  = 70000;
            nu = 0.3;

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

  

        function [LHSr,Mr, RHSr, lhs] = createElasticProblem(obj)
            u = LagrangianFunction.create(obj.meshDomain,obj.meshDomain.ndim,'P1');
            material = obj.createMaterial(obj.meshDomain);
            [lhs,LHSr] = obj.computeStiffnessMatrix(obj.meshDomain,u,material);
            [~,Mr] = obj.computeMassMatrix(obj.meshDomain,u,material);
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

        function [M,Mr] = computeMassMatrix(obj,mesh,dispFun,mat)
            s.type     = 'MassMatrix';
            s.mesh     = mesh;
            s.test      = dispFun;
            s.trial      = dispFun;
            s.material = mat;
            %s.quadratureOrder = 'LINEAR';
            lhs = LHSIntegrator.create(s);
            M   = lhs.compute();
            Mr = obj.bcApplier.fullToReducedMatrixDirichlet(M);
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
            EIFEMfilename = obj.fileNameEIFEM;
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
            d = DomainDecompositionDofManager3D(s);
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
