classdef DataGenerationPoissonFetiDPV2 < handle    
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
        function obj = DataGenerationPoissonFetiDPV2()
            close all;
            obj.numSubdomains = [15 15]; 
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
            obj.fetiSolver = FetiDPPoissonV2(obj.globalMesh, obj.localMeshes, obj.localStiffness, obj.localForces, obj.nodeTol, obj.dirichletDofs);
            
            % Visulatizació dels nodes
            obj.fetiSolver.visualizeFetiNodes();
            
            % 1. Comparativa de convergència global (3 casos)
            disp('--- Starting Global Convergence Comparison (3 Cases) ---');
            tic;
            [uMono, uFeti] = obj.runConvergenceComparison();
            tiempoComparativa = toc; 
            fprintf('Global comparison time: %.4f seconds\n\n', tiempoComparativa);

            % 2. Execució aïllada i recàlcul de FETI-DP + Dirichlet
             %obj.runIsolatedFetiDirichlet();

            % Visualització
            obj.visualizePoissonSolution(uFeti, 'FETI-DP Solution (Poisson)');
            obj.visualizePoissonSolution(uMono, 'Monolithic Solution (Poisson)');

            % Comparació d'error
            relError = norm(uFeti - uMono) / norm(uMono);
            fprintf('Relative error between FETI-DP and Monolithic: %e\n', relError);

            if relError < 1e-10
                disp('FETI-DP solution matches the direct solver.');
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
            disp('--- Assembling Global Monolithic System (Poisson) ---');
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
            fprintf('Maximum scalar value (Monolithic) = %.4e\n', maxUyMono);
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
            c.Label.String = 'Scalar Field u';
            title(titleStr);
            xlabel('X'); ylabel('Y');
        end
        function mS = createStructuredMesh(~, nPerSide)
            if nargin < 2
                nPerSide = 8;
            end
            x1 = linspace(0, 1/5, nPerSide);
            x2 = linspace(0, 1/5, nPerSide);
            [xv, yv] = meshgrid(x1, x2);
            [F, V] = mesh2tri(xv, yv, zeros(size(xv)), 'x');
            s.coord  = V(:,1:2);
            s.connec = F;
            s.interpType = 'LINEAR';
            mS = Mesh.create(s);
        end
       
        function [uMono, uFeti] = runConvergenceComparison(obj)
            tol = 1e-10;
        
            % --- CAS 1: CG i Càlcul Monolític ---
            uGlobal   = LagrangianFunction.create(obj.globalMesh, 1, 'P1');
            weakFormK = @(u,v) DP(Grad(v), Grad(u));
            kGlobal   = IntegrateLHS(weakFormK, uGlobal, uGlobal, obj.globalMesh, 'Domain', 2);
            weakFormM = @(u,v) DP(v, u);
            mGlobal   = IntegrateLHS(weakFormM, uGlobal, uGlobal, obj.globalMesh, 'Domain', 2);
            fGlobal   = mGlobal * ones(uGlobal.nDofs, 1);
            
            freeDofs  = setdiff(1:uGlobal.nDofs, obj.dirichletDofs);
            kFree     = kGlobal(freeDofs, freeDofs);
            fFree     = fGlobal(freeDofs);
            
            % Càlcul de la solució directa Monolítica
            uRed = kFree \ fFree;
            uMono = zeros(uGlobal.nDofs, 1);
            uMono(freeDofs) = uRed;
            
            x0mono    = zeros(length(fFree), 1);
            Pid       = @(r) r;

            tic;
            [~, residual1, ~, ~] = PCG.solve(@(x) kFree*x, fFree, x0mono, Pid, tol, uRed);
            timeMono = toc;

            kappaK = condest(kFree);
            nDofsMono = length(fFree);
        
            % --- CAS 2: CG Interfície FETI-DP ---
            [fMat, dBar] = obj.fetiSolver.assembleProblem();
            x0feti    = zeros(size(dBar));
            
            % Solució exacta del multiplicador (per fer anàlisis d'error al PCG si cal)
            lambdaExact = fMat \ dBar; 
            
            tic;
            [~, residual2, ~, ~] = PCG.solve(@(x) fMat*x, dBar, x0feti, Pid, tol, lambdaExact);
            timeFetiDual = toc;

            kappaF = condest(fMat);
            nDofsFeti = length(dBar);
        
            % --- CAS 3: PCG FETI-DP + Dirichlet ---
            Pdir = @(r) obj.fetiSolver.applyDirichletPrecond(r);
        
            % Recuperem la solució lambda del mètode iteratiu
            tic;
            [lambdaFetiPCG, residual3, ~, ~] = PCG.solve(@(x) fMat*x, dBar, x0feti, Pdir, 1e-10, lambdaExact);
            timeFetiDir = toc;

            M = obj.fetiSolver.buildPrecondMatrix();
            kappaPCG = condest(M * fMat);
        
            % Comparative Table
            obj.printComparisonTable(...
                length(residual1), length(residual2), length(residual3), ...
                kappaK, kappaF, kappaPCG, ...
                nDofsMono, nDofsFeti, ...
                prod(obj.numSubdomains), timeMono, timeFetiDual, timeFetiDir);
        
            % Graphics
            obj.plotConvergenceComparison(residual1, residual2, residual3, tol);
            
            % Reconstruction of uFeti using lambda solution from PCG
            uFeti = obj.fetiSolver.reconstructGlobalSolution(lambdaFetiPCG, obj.globalMesh.nnodes);
        end
        
        function printComparisonTable(~, iter1, iter2, iter3, kK, kF, kPCG, ...
                               nDofsMono, nDofsFeti, nSub, t1, t2, t3)
            fprintf('\n');
            fprintf('╔══════════════════════════════╦══════════╦══════════════╦══════════╦═════════════╗\n');
            fprintf('║ Case                         ║  Iter.   ║  κ (cond.)   ║  DOFs    ║ CPU Time(s) ║\n');
            fprintf('╠══════════════════════════════╬══════════╬══════════════╬══════════╬═════════════╣\n');
            fprintf('║ Monolithic CG (K·u = f)      ║  %5d     ║  %8.6e       ║  %6d     ║  %8.4f   ║\n', iter1, kK,   nDofsMono, t1);
            fprintf('║ FETI-DP Interface CG         ║  %5d     ║  %8.6e       ║  %6d     ║  %8.4f   ║\n', iter2, kF,   nDofsFeti, t2);
            fprintf('║ FETI-DP PCG + Dirichlet      ║  %5d     ║  %8.6e       ║  %6d     ║  %8.4f   ║\n', iter3, kPCG, nDofsFeti, t3);
            fprintf('╚══════════════════════════════╩══════════╩══════════════╩══════════╩═════════════╝\n');
            fprintf('  Subdomains: %d×%d = %d\n', round(sqrt(nSub)), round(sqrt(nSub)), nSub);
            fprintf('\n');
        end
        
        function plotConvergenceComparison(~, h1, h2, h3, tol)
            figure('Name', 'CG Convergence - Comparison', ...
                   'Color', 'w', 'Position', [100 100 750 480]);
        
            semilogy(1:length(h1), h1, '-o',  ...
                     'Color', [0.00 0.45 0.74], ...
                     'LineWidth', 1.8, 'MarkerSize', 4, ...
                     'DisplayName', 'Monolithic CG');
            hold on;
        
            semilogy(1:length(h2), h2, '-s',  ...
                     'Color', [0.85 0.33 0.10], ...
                     'LineWidth', 1.8, 'MarkerSize', 4, ...
                     'DisplayName', 'FETI-DP Interface CG');
        
            semilogy(1:length(h3), h3, '-^',  ...
                     'Color', [0.47 0.67 0.19], ...
                     'LineWidth', 1.8, 'MarkerSize', 4, ...
                     'DisplayName', 'FETI-DP PCG + Dirichlet Preconditioner');
        
            yline(tol, '--k', 'LineWidth', 1.2, ...
                  'Label', sprintf('tol = %.0e', tol), ...
                  'LabelHorizontalAlignment', 'left');
        
            xlabel('Iteration', 'FontSize', 12);
            ylabel('Relative residual', 'FontSize', 12);
            title('CG Convergence - Comparativa Global', 'FontSize', 13);
            legend('Location', 'northeast', 'FontSize', 11);
            grid on;
            ax = gca; ax.GridAlpha = 0.3; ax.YMinorGrid = 'on';
            xlim([1, max([length(h1), length(h2), length(h3)]) + 1]);
            ylim([tol * 0.1, 2]);
            hold off;
        end
       
    end
end