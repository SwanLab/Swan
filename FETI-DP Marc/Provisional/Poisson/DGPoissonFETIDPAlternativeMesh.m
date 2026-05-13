classdef DGPoissonFETIDPAlternativeMesh < handle    
    properties (Access = private)
        globalMesh       
        localMeshes      
        numSubdomains    
        nodeTol          
        
        dirichletDofs   
        localStiffness   
        localForces      
        fetiSolver       
    end
    
    methods (Access = public)
        function obj = DGPoissonFETIDPAlternativeMesh()
            close all;
            obj.numSubdomains = [5 5]; 
            obj.nodeTol = 1e-10;
            
            % Malla global i locals
            referenceMesh = obj.createStructuredMesh();
            s.nsubdomains   = obj.numSubdomains;
            s.meshReference = referenceMesh;
            s.tolSameNode   = obj.nodeTol;
            m = MeshCreatorFromRVE.create(s);
            [obj.globalMesh, obj.localMeshes, ~, ~, ~, ~, ~] = m.create();
            
            % Condicions de contorn de Dirichlet
            obj.computeDirichletDofs();
            
            % Càlcul de LHS i RHS
            obj.computeLHS();
            obj.computeRHS();
            
            % Solver FETI-DP
            obj.fetiSolver = FetiDPPoissonAM(obj.globalMesh, obj.localMeshes, obj.localStiffness, obj.localForces, obj.nodeTol, obj.dirichletDofs);
            
            % Visulatizació dels nodes
            obj.fetiSolver.visualizeFetiNodes();
            
            % Comparativa de convergència
            obj.runConvergenceComparison();
            %obj.runScalabilityStudy(); 
            %obj.runSubdomainScalabilityStudy();

            % Ensamblatge i resolució de FETI-DP
            [fMat, dBar] = obj.fetiSolver.assembleProblem();

            lambdaSol = fMat \ dBar;

            x0 = zeros(size(dBar, 1), 1);
            fMatFun = @(x) fMat * x;
            %P = @(r) r;
            P = @(r) obj.fetiSolver.applyDirichletPrecond(r);

            disp('Iniciant ensamblatge FETI-DP...');
            tic; 
            [lambdaSol2, residual, err,errAnorm] = PCG.solve(fMatFun, dBar, x0, P, 1e-8, lambdaSol);
            tiempoEnsamblaje = toc; 
            fprintf('Temps de ensamblatge FETI-DP: %.4f segons\n', tiempoEnsamblaje);


            % Reconstrucció
            uFeti = obj.fetiSolver.reconstructGlobalSolution(lambdaSol2, obj.globalMesh.nnodes);
            
            obj.visualizePoissonSolution(uFeti, 'Solució FETI-DP (Poisson)');
            maxUFeti = max(abs(uFeti));
            fprintf('Màxim valor u (FETI-DP) = %.4e\n', maxUFeti);
            
            % Comparació monolítica
            uMono = obj.solveMonolithicPoisson();
            obj.visualizePoissonSolution(uMono, 'Solució Monolítica (Poisson)');
            
            relError = norm(uFeti - uMono) / norm(uMono);
            fprintf('Error relatiu entre FETI-DP i Monolític: %e\n', relError);
            
            if relError < 1e-10
                disp('La solució FETI-DP coincideix amb el solver directe.');
            end           
        end
    end
    
    methods (Access = private)
        
        function computeDirichletDofs(obj)
            coords = obj.globalMesh.coord;

            minX = min(coords(:,1)); 
            maxX = max(coords(:,1));

            minY = min(coords(:,2)); 
            maxY = max(coords(:,2));
            
            isDirNode = (abs(coords(:,1) - minX) < obj.nodeTol) | ...
                        (abs(coords(:,1) - maxX) < obj.nodeTol) | ...
                        (abs(coords(:,2) - minY) < obj.nodeTol) | ...
                        (abs(coords(:,2) - maxY) < obj.nodeTol);
                        
            % En Poisson hi ha 1 DoF per Node, per tant coincideixen
            obj.dirichletDofs = find(isDirNode);
        end

        function computeLHS(obj)
            % Construeix únicament les matrius de rigidesa (K)
            numSub = prod(obj.numSubdomains);
            kCell = cell(numSub, 1);

            for i = 1:numSub
                localMesh = obj.localMeshes{i};
                uLoc = LagrangianFunction.create(localMesh, 1, 'P1');

                weakFormK = @(u,v) DP(Grad(v), Grad(u));
                kCell{i} = IntegrateLHS(weakFormK, uLoc, uLoc, localMesh, 'Domain', 2);
            end

            obj.localStiffness = kCell;
        end
        
        function computeRHS(obj)
            % Construeix únicament els vectors de força (F)
            numSub = prod(obj.numSubdomains);
            fCell = cell(numSub, 1);
            
            for i = 1:numSub
                localMesh = obj.localMeshes{i};
                uLoc = LagrangianFunction.create(localMesh, 1, 'P1');
                
                weakFormM = @(u,v) DP(v, u);
                mLocal = IntegrateLHS(weakFormM, uLoc, uLoc, localMesh, 'Domain', 2);
                fLocVec = mLocal * ones(uLoc.nDofs, 1);
                
                fCell{i} = fLocVec;
            end
            
            obj.localForces = fCell;
        end
        
        function uDirect = solveMonolithicPoisson(obj)
            disp('--- Ensamblant el Sistema Monolític Global (Poisson) ---');
            uGlobal = LagrangianFunction.create(obj.globalMesh, 1, 'P1');
            
            % LHS
            weakFormK = @(u,v) DP(Grad(v), Grad(u));
            kGlobal = IntegrateLHS(weakFormK, uGlobal, uGlobal, obj.globalMesh, 'Domain', 2);
            
            % RHS
            weakFormM = @(u,v) DP(v, u);
            mGlobal = IntegrateLHS(weakFormM, uGlobal, uGlobal, obj.globalMesh, 'Domain', 2);
            fGlobal = mGlobal * ones(uGlobal.nDofs, 1);
            
            % Apliquem els DoFs guardats a l'objecte, evitem recalcular
            freeDofs = setdiff(1:uGlobal.nDofs, obj.dirichletDofs);
            
            kRed = kGlobal(freeDofs, freeDofs);
            fRed = fGlobal(freeDofs);
            uRed = kRed \ fRed;
            
            uDirect = zeros(uGlobal.nDofs, 1);
            uDirect(freeDofs) = uRed;
            
            maxUyMono = max(abs(uDirect));
            fprintf('Màxim valor escalar (Monolític) = %.4e\n', maxUyMono);
        end
        
        function visualizePoissonSolution(obj, uGlobal, titleStr)
            coords = obj.globalMesh.coord;
            connec = obj.globalMesh.connec;
        
            figure('Name', titleStr, 'Color', 'w');
            hold on; axis equal;
        
            patch('Faces', connec, 'Vertices', coords, ...
                  'FaceVertexCData', uGlobal, 'FaceColor', 'interp', ...
                  'EdgeColor', '#333333', 'LineWidth', 0.5);
        
            colormap(jet);
            c = colorbar;
            c.Label.String = 'Camp Escalar u';
            title(titleStr);
            xlabel('X'); ylabel('Y');
        end
        
        function mS = createStructuredMesh(~)
            % Carrega la malla del fitxer .mat
            data  = load('DEF_Q4auxL_1.mat');
            e     = data.EIFEoper;
            mesh  = e.MESH;
        
            coord = mesh.COOR;           
            cnQ4  = double(mesh.CN);    
        
            % Conversion Q4 to triangles (2 triangles for quadrilateral element)
            tri1   = cnQ4(:, [1 2 3]);
            tri2   = cnQ4(:, [1 3 4]);
            connec = [tri1; tri2];       
        
            s.coord      = coord;
            s.connec     = connec;
            s.interpType = 'LINEAR';
            mS = Mesh.create(s);
        end

        % function mS = createStructuredMesh(~, nPerSide)
        %     if nargin < 2
        %         nPerSide = 8;
        %     end
        %     x1 = linspace(0, 1/5, nPerSide);
        %     x2 = linspace(0, 1/5, nPerSide);
        %     [xv, yv] = meshgrid(x1, x2);
        %     [F, V] = mesh2tri(xv, yv, zeros(size(xv)), 'x');
        % 
        %     s.coord  = V(:,1:2);
        %     s.connec = F;
        %     s.interpType = 'LINEAR';
        %     mS = Mesh.create(s);
        % end

        % function runConvergenceComparison(obj)
        %     tol     = 1e-10;
        % 
        %     % CAS 1: CG sobre sistema monolític global
        %     uGlobal   = LagrangianFunction.create(obj.globalMesh, 1, 'P1');
        % 
        %     weakFormK = @(u,v) DP(Grad(v), Grad(u));
        %     kGlobal   = IntegrateLHS(weakFormK, uGlobal, uGlobal, obj.globalMesh, 'Domain', 2);
        % 
        %     weakFormM = @(u,v) DP(v, u);
        %     mGlobal   = IntegrateLHS(weakFormM, uGlobal, uGlobal, obj.globalMesh, 'Domain', 2);
        %     fGlobal   = mGlobal * ones(uGlobal.nDofs, 1);
        % 
        %     freeDofs  = setdiff(1:uGlobal.nDofs, obj.dirichletDofs);
        %     kFree      = kGlobal(freeDofs, freeDofs);
        %     fFree     = fGlobal(freeDofs);
        % 
        %     x0mono    = zeros(length(fFree), 1);
        %     Afun1     = @(x) kFree * x;
        %     Pidentity       = @(r) r;
        %     xsol_zero = zeros(size(fFree)); % xsol desconeguda → zeros (err no es fa servir)
        % 
        %     [~, residual1, ~, ~] = PCG.solve(Afun1, fFree, x0mono, Pidentity, tol, xsol_zero);
        %     fprintf('CAS 1: CG Monolític → %d iteracions\n', length(residual1));
        % 
        %     % CAS 2: CG interfície FETI-DP (sense precondicionador)
        %     [fMat, dBar] = obj.fetiSolver.assembleProblem();
        %     x0feti    = zeros(size(dBar));
        %     Afun2     = @(x) fMat * x;
        %     xsol_feti = zeros(size(dBar));
        % 
        %     [~, residual2, ~, ~] = PCG.solve(Afun2, dBar, x0feti, Pidentity, tol, xsol_feti);
        %     fprintf('CAS 2: CG Interfície FETI → %d iteracions\n', length(residual2));
        % 
        %     % CAS 3: PCG interfície FETI-DP + Dirichlet preconditioner
        %     Pdir = @(r) obj.fetiSolver.applyDirichletPrecond(r, fMat);
        % 
        %     [~, residual3, ~, ~] = PCG.solve(Afun2, dBar, x0feti, Pdir, tol, xsol_feti);
        %     fprintf('CAS 3: PCG FETI + Dirichlet → %d iteracions\n', length(residual3));
        % 
        %     % Gràfica 
        %     obj.plotConvergenceComparison(residual1, residual2, residual3, tol);
        % end

        function runConvergenceComparison(obj)
            tol = 1e-10;
        
            % CASE 1: CG Monolithic
            uGlobal   = LagrangianFunction.create(obj.globalMesh, 1, 'P1');
            weakFormK = @(u,v) DP(Grad(v), Grad(u));
            kGlobal   = IntegrateLHS(weakFormK, uGlobal, uGlobal, obj.globalMesh, 'Domain', 2);
            weakFormM = @(u,v) DP(v, u);
            mGlobal   = IntegrateLHS(weakFormM, uGlobal, uGlobal, obj.globalMesh, 'Domain', 2);
            fGlobal   = mGlobal * ones(uGlobal.nDofs, 1);
        
            freeDofs  = setdiff(1:uGlobal.nDofs, obj.dirichletDofs);
            kFree     = kGlobal(freeDofs, freeDofs);
            fFree     = fGlobal(freeDofs);
        
            x0mono    = zeros(length(fFree), 1);
            Pid       = @(r) r;
            xsol_zero = zeros(size(fFree));
        
            [~, residual1, ~, ~] = PCG.solve(@(x) kFree*x, fFree, x0mono, Pid, tol, xsol_zero);
        
            % Condition Number case 1
            kappaK = condest(kFree);
            nDofsMono = length(fFree);
        
            % CASE 2: CG Interface FETI-DP 
            [fMat, dBar] = obj.fetiSolver.assembleProblem();
            x0feti    = zeros(size(dBar));
            xsol_feti = zeros(size(dBar));
        
            [~, residual2, ~, ~] = PCG.solve(@(x) fMat*x, dBar, x0feti, Pid, tol, xsol_feti);
        
            % Condition Number case 2
            kappaF = condest(fMat);
            nDofsFeti = length(dBar);
        
            % CASE 3: PCG FETI-DP + Dirichlet
            Pdir = @(r) obj.fetiSolver.applyDirichletPrecond(r);
        
            [~, residual3, ~, ~] = PCG.solve(@(x) fMat*x, dBar, x0feti, Pdir, tol, xsol_feti);
        
            % Condition Number case 3: κ(M^-1 F)
            M = obj.fetiSolver.buildPrecondMatrix();
            MInvF = M \ full(fMat);
            kappaPCG = condest(MInvF);
        
            % Table
            obj.printComparisonTable(...
                length(residual1), length(residual2), length(residual3), ...
                kappaK, kappaF, kappaPCG, ...
                nDofsMono, nDofsFeti, ...
                prod(obj.numSubdomains));
      
            obj.plotConvergenceComparison(residual1, residual2, residual3, tol);
        end
        
        function printComparisonTable(~, iter1, iter2, iter3, kK, kF, kPCG, ...
                                       nDofsMono, nDofsFeti, nSub)
            fprintf('\n');
            fprintf('╔════════════════╦══════════╦══════════════╦══════════╗\n');
            fprintf('║ Case                         ║  Iter.   ║  κ (cond.)   ║  DOFs    ║\n');
            fprintf('╠════════════════╬══════════╬══════════════╬══════════╣\n');
            fprintf('║ CG Monolithic  (K·u = f)     ║  %5d   ║  %8.6e       ║  %6d     ║\n', iter1, kK,   nDofsMono);
            fprintf('║ CG Interface FETI-DP         ║  %5d   ║  %8.6e       ║  %6d     ║\n', iter2, kF,   nDofsFeti);
            fprintf('║ PCG FETI-DP + Dirichlet      ║  %5d   ║  %8.6e       ║  %6d     ║\n', iter3, kPCG, nDofsFeti);
            fprintf('╚════════════════╩══════════╩══════════════╩══════════╝\n');
            fprintf('  Subdominis: %d×%d = %d\n', round(sqrt(nSub)), round(sqrt(nSub)), nSub);
            fprintf('\n');
        end

        function plotConvergenceComparison(~, h1, h2, h3, tol)
            figure('Name', 'Convergència CG - Comparativa', ...
                   'Color', 'w', 'Position', [100 100 750 480]);
        
            semilogy(1:length(h1), h1, '-o',  ...
                     'Color', [0.00 0.45 0.74], ...
                     'LineWidth', 1.8, 'MarkerSize', 4, ...
                     'DisplayName', 'CG Monolític');
            hold on;
        
            semilogy(1:length(h2), h2, '-s',  ...
                     'Color', [0.85 0.33 0.10], ...
                     'LineWidth', 1.8, 'MarkerSize', 4, ...
                     'DisplayName', 'CG Interfície FETI-DP');
        
            semilogy(1:length(h3), h3, '-^',  ...
                     'Color', [0.47 0.67 0.19], ...
                     'LineWidth', 1.8, 'MarkerSize', 4, ...
                     'DisplayName', 'PCG FETI-DP + Dirichlet Preconditioner');
        
            % Línia de tolerància
            yline(tol, '--k', 'LineWidth', 1.2, ...
                  'Label', sprintf('tol = %.0e', tol), ...
                  'LabelHorizontalAlignment', 'left');
        
            xlabel('Iteració', 'FontSize', 12);
            ylabel('Residu relatiu', 'FontSize', 12);
            %ylabel('Residu relatiu  ||r_k|| / ||r_0||', 'FontSize', 12);
            title('Convergència CG - Malla 10×10 subdominis', 'FontSize', 13);
            legend('Location', 'northeast', 'FontSize', 11);
            grid on;
            ax = gca;
            ax.GridAlpha = 0.3;
            ax.YMinorGrid = 'on';
            xlim([1, max([length(h1), length(h2), length(h3)]) + 1]);
            ylim([tol * 0.1, 2]);
            hold off;
        end

    end
end