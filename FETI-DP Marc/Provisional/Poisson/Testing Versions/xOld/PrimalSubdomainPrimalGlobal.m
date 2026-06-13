classdef PrimalSubdomainPrimalGlobal < handle

    properties (Access = public)
        primalDofsGlobal
        dualDofsGlobal    
        remDofsGlobal

        K_local  
        f_local
    end
    properties (Access = private)
        meshDomain
        boundaryConditions
        bcApplier
        LHS
        RHS

        fileNameEIFEM
        tolSameNode
        nSubdomains

        Bd_local

    end

    methods (Access = public)

        function obj = PrimalSubdomainPrimalGlobal()
            close all
            obj.init()

            mR = obj.createReferenceMesh();
            bS  = mR.createBoundaryMesh();
            [mD,mSb,iC,lG,iCR,discMesh] = obj.createMeshDomain(mR);
            obj.meshDomain = mD;

            obj.extractFETIdofs(mSb, lG);

            obj.visualizeFETINodes();

            obj.computeLocalMatrices(mSb);

            [F, dbar] = assembleF_and_d(obj);
            
            % 1. SOLUCIÓ DIRECTA (DE REFERENCIA)
            tic;
            lambda_direct = F \ dbar;
            toc;
            
            % 2. SOLUCIÓ ITERATIVA FETI-DP (PCG + PRECONDICIONADOR)
            tol = 1e-12;
            maxit = 1000;
            
            % Definimos la función anónima que MATLAB usará como Precondicionador
            M_dirichlet = @(r) obj.applyDirichletPrecond(r);
            
            tic
            % Llamada al solver iterativo
            [lambda_sol, flag, relres, iter, resvec] = pcg(F, dbar, tol, maxit, M_dirichlet);
            toc
            
            if flag == 0
                disp(['¡Convergencia alcanzada en ', num2str(iter), ' iteraciones!']);
            else
                disp(['Atención: El solver terminó con flag ', num2str(flag)]);
            end
            
            % Gráfica de convergencia
            figure;
            semilogy(0:iter, resvec/resvec(1), '-o', 'LineWidth', 2, 'MarkerFaceColor', 'b');
            title('Convergencia de FETI-DP (Precondicionador de Dirichlet)');
            xlabel('Iteraciones');
            ylabel('Residuo Relativo ||r|| / ||r_0||');
            grid on;
            set(gca, 'YScale', 'log');

            % --- 3. RECONSTRUCCIÓN DE LA SOLUCIÓN GLOBAL Y VISUALIZACIÓN ---
            disp('Reconstruyendo el vector de desplazamientos (U_full)...');
            tic;
            U_full = obj.reconstructGlobalSolution(lambda_sol);
            toc;
            
            % Factor de escala (si el desplazamiento es 0.001, con un factor de 10 se verá como 0.01)
            factor_escala = 5; 
            disp('Dibujando malla deformada...');
            obj.visualizeDeformedMesh(U_full, factor_escala);
            uyMax = max(abs(U_full(2:2:end)));
            fprintf('Max Uy (FETI-DP) = %.4e\n', uyMax);
            
            % ========================================================
            % 4. COMPARACIÓN CON SOLVER MONOLÍTICO ("A MUERTE")
            % ========================================================
            U_mono = obj.solveMonolithicDirect();
            
            % Calcular error relativo entre ambos vectores de desplazamiento
            error_relativo = norm(U_full - U_mono) / norm(U_mono);
            fprintf('Error relativo entre FETI-DP y Monolítico: %e\n', error_relativo);
            
            if error_relativo < 1e-10
                disp('¡ÉXITO! FETI-DP y el solver monolítico coinciden matemáticamente.');
            else
                disp('Hay discrepancias entre los métodos.');
            end           
            
  
%             [bC,dir] = obj.createBoundaryConditions(obj.meshDomain);
%             obj.boundaryConditions = bC;
%             obj.createBCapplier()
% 
%             tic
%             [LHS,RHS,LHSf] = obj.createElasticProblem();
%             toc
%             obj.LHS = LHSf;
%             %             LHS = 0.5*(LHS+LHS');
% 
%             LHSf = @(x) LHS*x;
%             RHSf = RHS;
%             rhs2 = repmat(RHS,[1,20,1]);
%             tic
%             Usol = LHS\RHS;
%             toc
% %             Ufull = obj.bcApplier.reducedToFullVectorDirichlet(Usol);
%             %obj.plotSolution(Ufull,obj.meshDomain,1,1,0,obj.bcApplier,0)
% 
%             RBbasisFree = obj.forAlgebraicMultigrid();
%             Mid          = @(r) r;
%             Meifem       = obj.createEIFEMPreconditioner(mR,dir,iC,lG,bS,iCR,discMesh);
% %             MeifemCont   = obj.createEIFEMPreconditionerContinuous(mR,dir,iC,lG,bS,iCR,discMesh,obj.LHS);
%             Milu         = obj.createILUpreconditioner(LHS);
%             MgaussSeidel = obj.createGaussSeidelpreconditioner(LHS);
%             MJacobi      = obj.createJacobipreconditioner(LHS);
%             Mmodal       = obj.createModalpreconditioner(LHS);
% %            MblockD      = obj.createBlockDiagonalpreconditioner(LHS);
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
%             Mmult = @(r) Preconditioner.multiplePrec(r,LHSf,Milu,Meifem,Milu);
% % %            [eigVALMA_min,eigVALMA_max] = obj.computeEigs(LHS,Mmult);
% % % %            eigMA = sort([diag(eigVALMA_min);diag(eigVALMA_max)]);
% % %             eigMA = [diag(eigVALMA_min);diag(eigVALMA_max)];
% % %            [eigVALLHS_min,eigVALLHS_max] = obj.computeEigs(LHS);
% % % %            eigLHS = sort([diag(eigVALLHS_min);diag(eigVALLHS_max)]);
% % %             eigLHS = [diag(eigVALLHS_min);diag(eigVALLHS_max)];
% % %             figure
% % %             plot((eigMA),'o','MarkerFaceColor', 'b')
% % %             hold on 
% % %             plot((eigLHS), 'o', 'MarkerFaceColor', 'r')
% % %             xlabel('Number')
% % %             ylabel('Eigenvalue')
% % %             legend({'EIFEM','Coefficient matrix',},'FontSize',12)
% %             set(gca, 'YScale', 'log')
% 
% %              Mmult = @(r) Preconditioner.multiplePrec(r,Mid,Meifem,Mid,LHSf,RHSf,obj.meshDomain,obj.bcApplier);
% %             zmult = Mmult(r);
% 
% %             zfull = obj.bcApplier.reducedToFullVectorDirichlet(zmult);
%             %obj.plotSolution(zfull,obj.meshDomaopenin,0,0,2,obj.bcApplier,0)
% 
% %             zeifem = Meifem(r);
% %             zfull = obj.bcApplier.reducedToFullVectorDirichlet(zeifem);
%             %obj.plotSolution(zfull,obj.meshDomain,0,0,1,obj.bcApplier,0)
%            % x0 = zmult;
%             tic
%             %           tau = @(r,A) 1;
%             [uPCG,residualPCG,errPCG,errAnormPCG] = PCG.solve(LHSf,RHSf,x0,Mmult,tol,Usol);
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

        function U_direct = solveMonolithicDirect(obj)
            % SOLVEMONOLITHICDIRECT Resuelve el problema global ensamblando
            % directamente K y F de toda la malla ("a muerte").
            
            disp('--- Ensamblando Sistema Monolítico Global ---');
            
            % 1. Crear espacio de funciones y material para el dominio global
            u_global = LagrangianFunction.create(obj.meshDomain, obj.meshDomain.ndim, 'P1');
            material = obj.createMaterial(obj.meshDomain);
            
            % 2. Ensamblar Matriz de Rigidez Global (K) igual que en los subdominios
            C = material;
            f_weak = @(u,v) DDP(SymGrad(v), DDP(C, SymGrad(u)));
            K_global = IntegrateLHS(f_weak, u_global, u_global, obj.meshDomain, 'Domain', 2);
            
            % 3. Ensamblar Vector de Fuerzas Global (F)
            total_dofs = obj.meshDomain.nnodes * obj.meshDomain.ndim;
            F_global   = zeros(total_dofs, 1);
            
            tol        = obj.tolSameNode;
            maxx       = max(obj.meshDomain.coord(:,1));
            load_value = -2; % Misma carga que en computeLocalMatrices
            
            % Aplicar la carga en la cara derecha
            for j = 1:obj.meshDomain.nnodes
                if abs(obj.meshDomain.coord(j,1) - maxx) < tol
                    y_dof = (j - 1) * obj.meshDomain.ndim + 2;
                    F_global(y_dof) = F_global(y_dof) + load_value;
                end
            end
            
            % 4. Aplicar Condiciones de Dirichlet (Pared izquierda fija)
            minx = min(obj.meshDomain.coord(:,1));
            isDirichletNode = abs(obj.meshDomain.coord(:,1) - minx) < tol;
            dirichletNodes  = find(isDirichletNode);
            
            dirichletDofs = [];
            for d = 1:obj.meshDomain.ndim
                dirichletDofs = [dirichletDofs; (dirichletNodes - 1) * obj.meshDomain.ndim + d];
            end
            dirichletDofs = sort(dirichletDofs);
            
            % Grados de libertad libres (donde resolveremos)
            freeDofs = setdiff(1:total_dofs, dirichletDofs);
            
            % 5. Reducir el sistema y Resolver (U = K \ F)
            disp('Resolviendo U = K \ F ...');
            K_red = K_global(freeDofs, freeDofs);
            F_red = F_global(freeDofs);
            
            tic;
            U_red = K_red \ F_red;
            toc;
            
            % 6. Reconstruir el vector completo insertando ceros en Dirichlet
            U_direct = zeros(total_dofs, 1);
            U_direct(freeDofs) = U_red;
            
            uyMax = max(abs(U_direct(2:2:end)));
            fprintf('Max Uy (Monolítico) = %.4e\n', uyMax);
        end

    end

    methods (Access = private)

        function init(obj)
            obj.nSubdomains  = [2 2]; %nx ny
%             obj.fileNameEIFEM = 'DEF_Q4auxL_1.mat';
%             obj.fileNameEIFEM = 'DEF_auxNew_2.mat';
            obj.fileNameEIFEM = 'DEF_Q4porL_1_raul.mat';
            obj.tolSameNode = 1e-10;

        end

        function [mD,mSb,iC,lG,iCR,discMesh] = createMeshDomain(obj,mR)
            s.nsubdomains   = obj.nSubdomains; %nx ny
            s.meshReference = mR;
            s.tolSameNode = obj.tolSameNode;
            m = MeshCreatorFromRVE.create(s);
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
            x1      = linspace(0,1,3);
            x2      = linspace(0,1,3);
            [xv,yv] = meshgrid(x1,x2);
            [F,V]   = mesh2tri(xv,yv,zeros(size(xv)),'x');
            s.coord  = V(:,1:2);
            s.connec = F;
            bgMesh = Mesh.create(s);

            s.coord    = s.coord(:,1:2);
            s.connec   = F;
            s.interpType = 'LINEAR';
            mS         = Mesh.create(s);
        end



        function levelSet = createLevelSetFunction(obj,bgMesh)
            sLS.type        = 'CircleInclusion';
            sLS.xCoorCenter = 0.5;
            sLS.yCoorCenter = 0.5;
            sLS.radius      = 0.2;
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
            filename = obj.fileNameEIFEM;
            load(filename);
            s.coord    = EIFEoper.MESH.COOR;
            isMin = s.coord==min(s.coord);
            isMax = s.coord==max(s.coord);

            s.connec   = EIFEoper.MESH.CN;
            s.interType = 'QUADRATIC';
            mS         = Mesh.create(s);
        end

        function mCoarse = createCoarseMesh(obj,mR)
            s.nsubdomains   = obj.nSubdomains; %nx ny
            s.meshReference = obj.createReferenceCoarseMesh(mR);
%             s.meshReference = obj.loadReferenceCoarseMesh(mR);
            s.tolSameNode   = obj.tolSameNode;
            mRVECoarse      = MeshCreatorFromRVE.create(s);
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

        function cMesh = loadReferenceCoarseMesh(obj, mR)
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
            E  = 100000;
            %nu = 1/3;
%             E  = 70000;
            nu = 0.3;
            Epstr  = E/(1-nu^2);
            nupstr = nu/(1-nu);
            young   = ConstantFunction.create(Epstr,mesh);
            poisson = ConstantFunction.create(nupstr,mesh);
%             young   = ConstantFunction.create(E,mesh);
%             poisson = ConstantFunction.create(nu,mesh);
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

                        Dir{2}.domain    = @(coor) isRight(coor) ;
                        Dir{2}.direction = [1,2];
                        Dir{2}.value     = 0;

            PL.domain    = @(coor) isTop(coor);
            PL.direction = [2];
            PL.value     = [-0.1];
%                         PL.domain    = @(coor) isRight(coor);
%                         PL.direction = [1];
%                         PL.value     = [0.1];
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
            EIFEMfilename = obj.fileNameEIFEM;
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

        function Meifem = createEIFEMPreconditionerContinuous(obj,mR,dir,iC,lG,bS,iCR,dMesh,LHS)
            % obj.EIFEMfilename = '/home/raul/Documents/Thesis/EIFEM/RAUL_rve_10_may_2024/EXAMPLE/EIFE_LIBRARY/DEF_Q4porL_2s_1.mat';
            EIFEMfilename = obj.fileNameEIFEM;
            % obj.EIFEMfilename = '/home/raul/Documents/Thesis/EIFEM/05_HEXAG2D/EIFE_LIBRARY/DEF_Q4auxL_1.mat';
            filename        = EIFEMfilename;
            s.RVE           = TrainedRVE(filename);
            s.mesh          = obj.createCoarseMesh(mR);
%            s.mesh          = obj.loadCoarseMesh(mR);
            s.DirCond       = dir;
            s.nSubdomains = obj.nSubdomains;
            s.meshRef = mR;
            s.meshDomain = obj.meshDomain;
           
            eifem           = EIFEM_trying_ideas(s);


            ss.ddDofManager = obj.createDomainDecompositionDofManager(iC,lG,bS,mR,iCR);
            ss.EIFEMsolver = eifem;
            ss.bcApplier = obj.bcApplier;
            ss.dMesh     = dMesh;
            ss.LHS  = LHS; 
            ss.type = 'EIFEMcont';
            ss.nSubdomains = obj.nSubdomains;
            eP = Preconditioner.create(ss);
            k  = eP.computeKEIFEMglobal(LHS);
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
            refPoint = (min(obj.meshDomain.coord)+min(obj.meshDomain.coord))/2 ;

            RB = RigidBodyFunction.create(obj.meshDomain,refPoint);
            xt  = RB.basisFunctions{1}.project('P1');
            yt  = RB.basisFunctions{2}.project('P1');
            rot  = RB.basisFunctions{3}.project('P1');
            BG = [reshape(xt.fValues',[],1),reshape(yt.fValues',[],1),reshape(rot.fValues',[],1)];
            B = [obj.bcApplier.fullToReducedVectorDirichlet(BG(:,1)),obj.bcApplier.fullToReducedVectorDirichlet(BG(:,2)),...
                 obj.bcApplier.fullToReducedVectorDirichlet(BG(:,3))]; 
        end

        function [eigValsMin,eigValsMax ] = computeEigs(obj,A,M)
             opts = struct();
                opts.tol = 1e-7;          % Reasonably tight tolerance
                opts.maxit = 20000;        % More iterations than default
                opts.issym = true;        % Matrix is symmetric
                opts.isreal = true;       % Matrix is real
            if nargin == 2

                % Compute 6 largest-magnitude eigenvalues
                [eigVecs, eigValsMin] = eigs(A, 6809, 'smallestreal', opts);
                 [eigVecs, eigValsMax] = eigs(A, 1, 'largestreal', opts);

                % Extract eigenvalues from diagonal
%                 eigenvalues = diag(eigVals);
            else
                for j = 1:size(A,1)
                    MA(:, j) = M(A(:, j));
                end
                try
                    chol(MA);  % Should succeed if truly SPD
                    disp('Matrix is numerically SPD.');
                catch
                    disp('Matrix is not SPD numerically.');
                end
                [eigVecs, eigValsMin] = eigs(MA, 6809, 'smallestreal', opts);
                eigValsMin = real(eigValsMin);
                [eigVecs, eigValsMax] = eigs(MA, 1, 'largestreal', opts);
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

        % function extractFETIdofs(obj, mSb, ~)
        %     nSub = prod(obj.nSubdomains);
        %     pGlobal = cell(nSub, 1);
        %     dGlobal = cell(nSub, 1);
        %     rGlobal = cell(nSub, 1);
        % 
        %     tol = obj.tolSameNode;
        %     ndim = obj.meshDomain.ndim;
        %     globalCoords = obj.meshDomain.coord; 
        % 
        %     for i = 1:nSub
        %         coords = mSb{i}.coord;
        % 
        %         % Límites locales del subdominio
        %         minX = min(coords(:,1)); maxX = max(coords(:,1));
        %         minY = min(coords(:,2)); maxY = max(coords(:,2));
        % 
        %         % 1. Nodos Primales (Esquinas del subdominio local)
        %         isBL = abs(coords(:,1) - minX) < tol & abs(coords(:,2) - minY) < tol; % BottomLeft
        %         isBR = abs(coords(:,1) - maxX) < tol & abs(coords(:,2) - minY) < tol; % BottomRight
        %         isTL = abs(coords(:,1) - minX) < tol & abs(coords(:,2) - maxY) < tol; % TopLeft
        %         isTR = abs(coords(:,1) - maxX) < tol & abs(coords(:,2) - maxY) < tol; % TopRight
        % 
        %         localPrimalNodes = find(isBL | isBR | isTL | isTR);
        % 
        %         % 2. Identificar TODOS los nodos en la frontera del subdominio
        %         isOnLocalBoundary = abs(coords(:,1) - minX) < tol | abs(coords(:,1) - maxX) < tol | ...
        %                             abs(coords(:,2) - minY) < tol | abs(coords(:,2) - maxY) < tol;
        %         localBoundaryNodes = find(isOnLocalBoundary);
        % 
        %         % 3. Nodos Duales: Toda la frontera local (interiores y exteriores) MENOS los primales
        %         localDualNodes = setdiff(localBoundaryNodes, localPrimalNodes);
        % 
        %         % 4. Nodos Restantes (Remaining): Estrictamente los nodos internos al subdominio
        %         allLocalNodes = (1:size(coords,1))';
        %         localRemNodes = setdiff(allLocalNodes, localBoundaryNodes);
        % 
        %         % 5. Mapeo a índices globales y conversión a Grados de Libertad (DoFs)
        %         [~, globalNodes] = ismembertol(coords, globalCoords, tol, 'ByRows', true);                
        % 
        %         pGlobal{i} = obj.nodesToDofs(globalNodes(localPrimalNodes), ndim);
        %         dGlobal{i} = obj.nodesToDofs(globalNodes(localDualNodes), ndim);
        %         rGlobal{i} = obj.nodesToDofs(globalNodes(localRemNodes), ndim);
        %     end
        % 
        %     obj.primalDofsGlobal = pGlobal;
        %     obj.dualDofsGlobal   = dGlobal;
        %     obj.remDofsGlobal    = rGlobal;          
        % end
        
        function extractFETIdofs(obj, mSb, ~)
            nSub = prod(obj.nSubdomains);
            pGlobal = cell(nSub, 1);
            dGlobal = cell(nSub, 1);
            rGlobal = cell(nSub, 1);

            tol = obj.tolSameNode;
            ndim = obj.meshDomain.ndim;
            globalCoords = obj.meshDomain.coord; 

            % --- LIMITES GLOBALES DEL DOMINIO ---
            gMinX = min(globalCoords(:,1)); gMaxX = max(globalCoords(:,1));
            gMinY = min(globalCoords(:,2)); gMaxY = max(globalCoords(:,2));

            for i = 1:nSub
                coords = mSb{i}.coord;

                % Límites locales del subdominio
                minX = min(coords(:,1)); maxX = max(coords(:,1));
                minY = min(coords(:,2)); maxY = max(coords(:,2));

                % 1. Nodos Primales (Esquinas del subdominio local)
                isBL = abs(coords(:,1) - minX) < tol & abs(coords(:,2) - minY) < tol; % BottomLeft
                isBR = abs(coords(:,1) - maxX) < tol & abs(coords(:,2) - minY) < tol; % BottomRight
                isTL = abs(coords(:,1) - minX) < tol & abs(coords(:,2) - maxY) < tol; % TopLeft
                isTR = abs(coords(:,1) - maxX) < tol & abs(coords(:,2) - maxY) < tol; % TopRight

                localPrimalNodes = find(isBL | isBR | isTL | isTR);

                % 2. Identificar TODOS los nodos en la frontera del subdominio
                isOnLocalBoundary = abs(coords(:,1) - minX) < tol | abs(coords(:,1) - maxX) < tol | ...
                                    abs(coords(:,2) - minY) < tol | abs(coords(:,2) - maxY) < tol;

                % 3. Identificar fronteras GLOBALES LIBRES (Top, Bottom, Right)
                % Estas zonas NO deben ser Duales para que la malla pueda deformarse ahí
                isOnFreeGlobalBoundary = (abs(coords(:,1) - gMaxX) < tol) | ...
                                         (abs(coords(:,2) - gMinY) < tol) | ...
                                         (abs(coords(:,2) - gMaxY) < tol);

                % 4. Nodos Duales: Frontera local MENOS las fronteras libres globales
                % Esto deja únicamente las interfaces internas y la pared izquierda (que sí queremos fijar)
                isDual = isOnLocalBoundary & ~isOnFreeGlobalBoundary;
                localDualNodes = setdiff(find(isDual), localPrimalNodes);

                % 5. Nodos Restantes: Todo lo que no sea Dual (incluye interiores y fronteras libres)
                isRem = ~isDual;
                localRemNodes = setdiff(find(isRem), localPrimalNodes);

                % 6. Mapeo a índices globales y conversión a Grados de Libertad (DoFs)
                [~, globalNodes] = ismembertol(coords, globalCoords, tol, 'ByRows', true);                

                pGlobal{i} = obj.nodesToDofs(globalNodes(localPrimalNodes), ndim);
                dGlobal{i} = obj.nodesToDofs(globalNodes(localDualNodes), ndim);
                rGlobal{i} = obj.nodesToDofs(globalNodes(localRemNodes), ndim);
            end

            obj.primalDofsGlobal = pGlobal;
            obj.dualDofsGlobal   = dGlobal;
            obj.remDofsGlobal    = rGlobal;          
        end

        function dofs = nodesToDofs(~, nodes, ndim)            
            nodes = nodes(:);
            dofs = zeros(length(nodes)*ndim, 1);
            for d = 1:ndim
                dofs(d:ndim:end) = (nodes-1)*ndim + d;
            end
            dofs = sort(dofs);
        end

        function computeLocalMatrices(obj, mSb)
            nSub = prod(obj.nSubdomains);
            K_cell = cell(nSub, 1);
            f_cell = cell(nSub, 1);
            
            % Màxima coordenanda global per trobar la pared de la dreta
            globalCoords = obj.meshDomain.coord;
            maxx = max(globalCoords(:,1));
            tol = obj.tolSameNode;
            load_value = -2; % Càrrega

            for i = 1:nSub
                % Extreu local mesh
                localMesh = mSb{i};
                               
                % Crear espai de funcions local
                u_loc = LagrangianFunction.create(localMesh, localMesh.ndim, 'P1');
                
                % Poisson (para resolver poisson en lugar de elasticidad)
                %u_loc = LagrangianFunction.create(localMesh, 1, 'P1');
                %f_weak = @(u,v) DP(Grad(v),Grad(u));

                % Crear material
                mat_loc = obj.createMaterial(localMesh);
                
                % Ensamblar Matriu de Rigidessa Local K
                C = mat_loc;
                f_weak = @(u,v) DDP(SymGrad(v), DDP(C, SymGrad(u)));
                K_cell{i} = IntegrateLHS(f_weak, u_loc, u_loc, localMesh, 'Domain', 2);

                % Ensamblar Vector de Forces Local f
                f_loc_vec = zeros(u_loc.nDofs, 1);
            
                % Aplicar càrrega a la cara de la dreta
                localCoords = localMesh.coord;
                for j = 1:size(localCoords, 1)
                    if abs(localCoords(j,1) - maxx) < tol
                        
                        y_dof = (j - 1) * localMesh.ndim + 2; 
                        
                        % Aplicamos la carga hacia abajo
                        f_loc_vec(y_dof) = f_loc_vec(y_dof) + load_value;
                    end
                end
                
                f_cell{i} = f_loc_vec;
            end
            
            obj.K_local = K_cell;
            obj.f_local = f_cell;
        end

        function [Krr, Krp, Kpr, Kpp, fr, fp, Tdr] = splitLocalMatrices(obj, sub_id)
            K = obj.K_local{sub_id};
            f = obj.f_local{sub_id};
            
            % 1. Recuperamos los DoFs GLOBALES de este subdominio
            p_dofs_global = obj.primalDofsGlobal{sub_id};
            d_dofs_global = obj.dualDofsGlobal{sub_id};
            r_dofs_global = obj.remDofsGlobal{sub_id};
            
            % 2. Juntamos todos los globales que pertenecen a este subdominio
            all_dofs_global = sort([p_dofs_global; d_dofs_global; r_dofs_global]);
            
            % 3. Buscamos qué posición local (1 al N) ocupa cada DoF dentro de la matriz
            [~, p_dofs_local] = ismember(p_dofs_global, all_dofs_global);
            [~, d_dofs_local] = ismember(d_dofs_global, all_dofs_global);
            [~, r_dofs_local] = ismember(r_dofs_global, all_dofs_global);
            
            % 4. Agrupamos los "Restantes" de FETI (Internos + Duales)
            rem_dofs_local = sort([r_dofs_local; d_dofs_local]);
            p_dofs_local = sort(p_dofs_local);
            
            % 5. Particionamos la matriz K (Slicing)
            Krr = K(rem_dofs_local, rem_dofs_local);
            Krp = K(rem_dofs_local, p_dofs_local);
            Kpr = K(p_dofs_local, rem_dofs_local); % Igual a Krp'
            Kpp = K(p_dofs_local, p_dofs_local);
            
            % 6. Particionamos el vector f (Slicing)
            fr = f(rem_dofs_local);
            fp = f(p_dofs_local);
                        
            
            d_dofs_local_sorted = sort(d_dofs_local);
            
            nD = length(d_dofs_local_sorted); 
            nR = length(rem_dofs_local);      
            
            [~, pos_in_rem] = ismember(d_dofs_local_sorted, rem_dofs_local);
            
            rows = (1:nD)';
            cols = pos_in_rem;
            vals = ones(nD, 1);
            
            Tdr = sparse(rows, cols, vals, nD, nR);
        end

        function [F, dbar] = assembleF_and_d(obj)
            allPrimals = unique(vertcat(obj.primalDofsGlobal{:}));
            allDuals = setdiff(unique(vertcat(obj.dualDofsGlobal{:})), allPrimals);
            
            nP = length(allPrimals);
            nL = length(allDuals);
            
            Fdual = sparse(nL, nL);
            dbar   = zeros(nL, 1);
            
            SPP = sparse(nP, nP);
            BR_KRRinv_KRP = sparse(nL, nP);
            rhsPrimal = zeros(nP, 1);
            
            nSub = prod(obj.nSubdomains);
            
            visitedDuals = zeros(max(allDuals), 1); 
            
            for s_id = 1:nSub
                [Krr, Krp, Kpr, Kpp, fr, fp, Tdr] = obj.splitLocalMatrices(s_id);
                
                Urd = Krr \ full(Tdr');      % Krr^-1 * Tdr^T
                Urp = Krr \ full(Krp);       % Krr^-1 * Krp        
                Urf = Krr \ fr;             % Krr^-1 * f_r

                Ap = obj.createAp(s_id, allPrimals);
                [Bd, visitedDuals] = obj.createBd(s_id, allDuals, visitedDuals);
                obj.Bd_local{s_id} = Bd;
                                
                % ENSAMBLATGEE GLOBAL
                Spp_loc = Kpp - Kpr * Urp;
                SPP = SPP + Ap * Spp_loc * Ap';
                
                if size(Bd, 2) > 0
                    Fdual = Fdual + Bd * (Urd' * Tdr') * Bd';
                    dbar = dbar + Bd * (Urd' *  fr);
                    BR_KRRinv_KRP = BR_KRRinv_KRP + Bd * (Tdr * Urp) * Ap';
                end
                
                rhsPrimal = rhsPrimal + Ap * (fp - Urp' * fr);
            end

            % APLICAR CONDICIONS DE CONTORN 
            active = obj.getActivePrimalDofs(allPrimals);   
            
            SPP_a  = SPP(active, active);                    
            BR_a   = BR_KRRinv_KRP(:, active);               
            rhs_a  = rhsPrimal(active);                     
            
            F    = Fdual + BR_a * (SPP_a \ full(BR_a'));    
            dbar = dbar   - BR_a * (SPP_a \ rhs_a);         

            
            %F = F_dual + BR_KRRinv_KRP * (SPP \ full(BR_KRRinv_KRP'));
            %dbar = dbar - BR_KRRinv_KRP * (SPP \ rhs_primal);           
        end

        function Ap = createAp(obj, s_id, all_primals)
            p_global = obj.primalDofsGlobal{s_id};
            nP_global = length(all_primals);
            nP_local = length(p_global);
            
            [~, rows] = ismember(p_global, all_primals);
            cols = (1:nP_local)';
            vals = ones(nP_local, 1);
            
            Ap = sparse(rows, cols, vals, nP_global, nP_local);
        end
        
       function [Bd, visited_duals] = createBd(obj, s_id, all_duals, visited_duals)
            d_global = obj.dualDofsGlobal{s_id};
            nL_global = length(all_duals);
            nL_local  = length(d_global);
        
            if nL_local == 0
                Bd = sparse(nL_global, 0);
                return;
            end
        
            [~, rows] = ismember(d_global, all_duals);
            cols = (1:nL_local)';
        
            % +1 si és la primera visita, -1 si és la segona
            isFirstVisit   = (visited_duals(d_global) == 0);
            vals           = ones(nL_local, 1);
            vals(~isFirstVisit) = -1;
        
            % Actualitzem les visites dels dofs que eren 0
            visited_duals(d_global(isFirstVisit)) = 1;
        
            Bd = sparse(rows, cols, vals, nL_global, nL_local);
        end

        function active_idx = getActivePrimalDofs(obj, all_primals)
            coords = obj.meshDomain.coord;
            tol    = obj.tolSameNode;
            minx   = min(coords(:,1)); 
            
            nP = length(all_primals);
            isDirichlet = false(nP, 1);
            
            for k = 1:nP
                dof  = all_primals(k);
                node = ceil(dof / obj.meshDomain.ndim);
                x    = coords(node, 1);
                
                % Fixem pared esquerra
                if abs(x - minx) < tol
                    isDirichlet(k) = true;
                end
            end
            
            active_idx = find(~isDirichlet);  
        end

        function z = applyDirichletPrecond(obj, r)
            nL = length(r);
            z = zeros(nL, 1);
            nSub = prod(obj.nSubdomains);
            
            % Escalado simple (asumiendo multiplicidad 2 en interfaces internas)
            scaling = 0.5; 

            for s_id = 1:nSub
                Bd = obj.Bd_local{s_id};
                if size(Bd, 2) == 0
                    continue; % Si no tiene nodos duales, saltamos
                end
                
                K = obj.K_local{s_id};
                
                % Identificamos nodos duales (frontera) e internos (resto)
                d_global = obj.dualDofsGlobal{s_id};
                i_global = obj.remDofsGlobal{s_id}; 
                
                % Buscamos su posición en la matriz local
                all_dofs = sort([obj.primalDofsGlobal{s_id}; d_global; i_global]);
                [~, d_loc] = ismember(d_global, all_dofs);
                [~, i_loc] = ismember(i_global, all_dofs);
                
                % Extraemos los bloques de la matriz de rigidez
                K_ii = K(i_loc, i_loc);
                K_id = K(i_loc, d_loc);
                K_di = K(d_loc, i_loc);
                K_dd = K(d_loc, d_loc);
                
                % 1. Calcular el Complemento de Schur local (S_dd)
                S_local = K_dd - K_di * (K_ii \ K_id);
                
                % 2. Extraer el residuo local y escalarlo
                r_local = scaling * (Bd' * r);
                
                % 3. Resolver localmente (Schur * residuo)
                z_local = S_local * r_local;
                
                % 4. Ensamblar la contribución al vector global
                z = z + Bd * (scaling * z_local);
            end
        end


        function visualizeFETINodes(obj)
            % 1. Get global coordinates and number of dimensions
            globalCoords = obj.meshDomain.coord;
            ndim = obj.meshDomain.ndim;
            
            % 2. Gather all DOFs globally
            all_primal_dofs = unique(vertcat(obj.primalDofsGlobal{:}));
            all_dual_dofs   = unique(vertcat(obj.dualDofsGlobal{:}));
            all_rem_dofs    = unique(vertcat(obj.remDofsGlobal{:}));
            
            % 3. Convert DOFs back to Node IDs (Node = ceil(DoF / ndim))
            primal_nodes = unique(ceil(all_primal_dofs / ndim));
            dual_nodes   = unique(ceil(all_dual_dofs / ndim));
            rem_nodes    = unique(ceil(all_rem_dofs / ndim));
            
            % Remove overlaps to ensure clean visualization (Primal > Dual > Rem)
            dual_nodes = setdiff(dual_nodes, primal_nodes);
            rem_nodes  = setdiff(rem_nodes, union(primal_nodes, dual_nodes));
            
            % 4. Get the spatial coordinates for each set of nodes
            p_coords = globalCoords(primal_nodes, :);
            d_coords = globalCoords(dual_nodes, :);
            r_coords = globalCoords(rem_nodes, :);
            
            % 5. Plotting
            figure('Name', 'FETI-DP Nodes', 'Color', 'w');
            hold on; axis equal;
            
            % Plot the mesh wireframe in the background (if connectivity exists)
            if isprop(obj.meshDomain, 'connec')
                trimesh(obj.meshDomain.connec, globalCoords(:,1), globalCoords(:,2), ...
                    zeros(size(globalCoords(:,1))), 'EdgeColor', [0.8 0.8 0.8], 'FaceColor', 'none');
            end
            
            % Scatter plot for Remaining (Internal) nodes
            if ~isempty(r_coords)
                scatter(r_coords(:,1), r_coords(:,2), 20, [0.5 0.5 0.5], 'filled', 'DisplayName', 'Remaining (Internal)');
            end
            
            % Scatter plot for Dual (Interface) nodes
            if ~isempty(d_coords)
                scatter(d_coords(:,1), d_coords(:,2), 40, 'b', 'filled', 'DisplayName', 'Dual (Interface)');
            end
            
            % Scatter plot for Primal (Corner) nodes
            if ~isempty(p_coords)
                scatter(p_coords(:,1), p_coords(:,2), 80, 'r', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1, 'DisplayName', 'Primal (Corners)');
            end
            
            % Formatting
            legend('Location', 'bestoutside');
            title('FETI-DP Node Distribution');
            xlabel('X'); ylabel('Y');
            grid on;
            hold off;
        end

        function U_full = reconstructGlobalSolution(obj, lambda_sol)
            % RECONSTRUCTGLOBALSOLUTION Reconstruye el vector global de desplazamientos
            % a partir de los multiplicadores de Lagrange (lambda).
            
            % --- PARTE 1: Reconstrucción Primal (uP) ---
            allPrimals = unique(vertcat(obj.primalDofsGlobal{:}));
            nP = length(allPrimals);
            SPP = sparse(nP, nP);
            rhsPrimal = zeros(nP, 1);
            nSub = prod(obj.nSubdomains);
            
            % Guardamos estas matrices temporalmente para no recalcularlas en el paso 2
            Urp_cell = cell(nSub, 1);
            Urd_cell = cell(nSub, 1);
            Krr_inv_fr_cell = cell(nSub, 1);
            
            for s_id = 1:nSub
                [Krr, Krp, Kpr, Kpp, fr, fp, Tdr] = obj.splitLocalMatrices(s_id);
                Ap = obj.createAp(s_id, allPrimals);
                Bd = obj.Bd_local{s_id};
                
                % Resoluciones locales
                Urp = Krr \ full(Krp); 
                Urd = Krr \ full(Tdr');
                Krr_inv_fr = Krr \ fr;
                
                % Guardamos para la Parte 2
                Urp_cell{s_id} = Urp;
                Urd_cell{s_id} = Urd;
                Krr_inv_fr_cell{s_id} = Krr_inv_fr;
                
                % Ensamblaje del Schur primal global (SPP)
                Spp_loc = Kpp - Kpr * Urp;
                SPP = SPP + Ap * Spp_loc * Ap';
                
                % Término independiente primal modificado por lambda
                % rhs_loc = f_p + Urp' * (Tdr' * Bd' * lambda - f_r)
                term = Tdr' * (Bd' * lambda_sol) - fr;
                rhsPrimal = rhsPrimal + Ap * (fp + Urp' * term);
            end
            
            % Aplicar condiciones de Dirichlet al sistema primal global
            active = obj.getActivePrimalDofs(allPrimals);
            SPP_a = SPP(active, active);
            rhs_a = rhsPrimal(active);
            
            uP = zeros(nP, 1);
            uP(active) = SPP_a \ rhs_a; % Resolvemos el problema primal global
            
            % --- PARTE 2: Reconstrucción de los DoFs locales y mapeo a U_full ---
            ndim = obj.meshDomain.ndim;
            nnodes = obj.meshDomain.nnodes;
            U_full = zeros(nnodes * ndim, 1);
            
            for s_id = 1:nSub
                p_global = obj.primalDofsGlobal{s_id};
                d_global = obj.dualDofsGlobal{s_id};
                r_global = obj.remDofsGlobal{s_id};
                
                Ap = obj.createAp(s_id, allPrimals);
                Bd = obj.Bd_local{s_id};
                
                Urp = Urp_cell{s_id};
                Urd = Urd_cell{s_id};
                Krr_inv_fr = Krr_inv_fr_cell{s_id};
                
                % 1. Extraemos la parte primal que le toca a este subdominio
                uP_loc = Ap' * uP;
                
                % 2. Calculamos los desplazamientos restantes (internos + duales)
                % u_r = Krr^-1 * f_r - Urp * u_p - Urd * Bd' * lambda
                u_rem_loc = Krr_inv_fr - Urp * uP_loc - Urd * (Bd' * lambda_sol);
                
                % 3. Mapeamos de vuelta al vector global U_full
                U_full(p_global) = uP_loc;
                
                % En splitLocalMatrices, los dofs restantes se ordenan agrupando r y d
                all_rem_global_sorted = sort([r_global; d_global]);
                U_full(all_rem_global_sorted) = u_rem_loc;
            end
        end

        function visualizeDeformedMesh(obj, U_global, scaleFactor)
            % VISUALIZEDEFORMEDMESH Muestra la malla original y la deformada.
            % 
            % Inputs:
            %   U_global    - Vector columna con los desplazamientos globales.
            %                 Debe tener tamaño [nnodes * ndim, 1] y el formato
            %                 [ux1; uy1; ux2; uy2; ...].
            %   scaleFactor - (Opcional) Factor para exagerar la deformación.
            %                 Por defecto es 1.0 (escala real).
        
            if nargin < 3
                scaleFactor = 1.0;
            end
        
            % 1. Extraer coordenadas originales y conectividad de la malla
            coords = obj.meshDomain.coord;
            connec = obj.meshDomain.connec;
            ndim   = obj.meshDomain.ndim;
            nnodes = size(coords, 1);
        
            % Verificar que el tamaño de U_global es correcto
            if length(U_global) ~= nnodes * ndim
                error('El tamaño de U_global (%d) no coincide con nnodes * ndim (%d).', length(U_global), nnodes * ndim);
            end
        
            % 2. Reestructurar el vector 1D a una matriz 2D [nnodes, ndim]
            % Al usar reshape con ndim filas y transponer, mapeamos [ux; uy] correctamente
            U_reshaped = reshape(U_global, ndim, nnodes)';
        
            % 3. Calcular coordenadas deformadas
            def_coords = coords + scaleFactor * U_reshaped;
        
            % 4. Calcular la magnitud del desplazamiento para colorear la malla
            disp_mag = sqrt(U_reshaped(:,1).^2 + U_reshaped(:,2).^2);
        
            % 5. Crear la figura
            figure('Name', 'Malla Deformada FETI-DP', 'Color', 'w');
            hold on; axis equal;
        
            % Dibujar la malla original en gris claro (fondo)
            patch('Faces', connec, 'Vertices', coords, ...
                  'FaceColor', 'none', 'EdgeColor', [0.8 0.8 0.8], ...
                  'LineStyle', '--', 'DisplayName', 'Malla Original');
        
            % Dibujar la malla deformada coloreada por la magnitud del desplazamiento
            patch('Faces', connec, 'Vertices', def_coords, ...
                  'FaceVertexCData', disp_mag, 'FaceColor', 'interp', ...
                  'EdgeColor', '#333333', 'LineWidth', 0.5, 'DisplayName', 'Malla Deformada');
        
            % Configuración visual de la gráfica
            colormap(jet);
            c = colorbar;
            c.Label.String = 'Magnitud del Desplazamiento ||u||';
            
            title(sprintf('Solución FETI-DP (Factor de Escala: %g)', scaleFactor));
            xlabel('X'); ylabel('Y');
            legend('Location', 'bestoutside');
            grid on;
            hold off;
        end


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
