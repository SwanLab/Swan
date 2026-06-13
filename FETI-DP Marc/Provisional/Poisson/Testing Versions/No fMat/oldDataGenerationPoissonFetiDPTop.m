classdef DataGenerationPoissonFetiDPTop < handle    
    
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
        
        function obj = DataGenerationPoissonFetiDPTop()
            close all;
            
            % 1. Initialization Parameters
            obj.numSubdomains = [40 40]; 
            obj.nodeTol       = 1e-10;
            
            % 2. Global and Local Meshes Setup
            referenceMesh   = obj.createStructuredMesh();
            s.nsubdomains   = obj.numSubdomains;
            s.meshReference = referenceMesh;
            s.tolSameNode   = obj.nodeTol;
            m = MeshCreatorFromRVE.create(s);
            
            [obj.globalMesh, obj.localMeshes, ~, ~, ~, ~, ~] = m.create();
            
            % 3. Apply Boundary Conditions & Assemble System
            obj.computeDirichletDofs();
            obj.computeLHS();
            obj.computeRHS();
            
            % 4. Initialize FETI-DP Solver
            obj.fetiSolver = FetiDPPoissonTop(obj.globalMesh, obj.localMeshes, ...
                                             obj.localStiffness, obj.localForces, ...
                                             obj.nodeTol, obj.dirichletDofs);
            
            % 5. Global Convergence Comparison (3 Cases)
            disp('--- Starting Global Convergence Comparison (3 Cases) ---');
            tic;
            [uMono, uFeti] = obj.runConvergenceComparison();
            totalTime = toc; 
            fprintf('Global comparison time: %.4f seconds\n\n', totalTime);
            
            % 6. Error Analysis
            relError = norm(uFeti - uMono) / norm(uMono);
            fprintf('Relative error between FETI-DP and Monolithic: %e\n', relError);
            if relError < 1e-10
                disp('Success: FETI-DP solution matches the monolithic direct solver.');
            end          
        end
        
    end

    % =========================================================
    % PRIVATE METHODS: SETUP & INITIALIZATION
    % =========================================================
    methods (Access = private)
        
        function mS = createStructuredMesh(obj, nPerSide)
            if nargin < 2
                nPerSide = 12;
            end
            numSub = obj.numSubdomains(1);
            x1 = linspace(0, 1/numSub, nPerSide);
            x2 = linspace(0, 1/numSub, nPerSide);
            [xv, yv] = meshgrid(x1, x2);
            [F, V]   = mesh2tri(xv, yv, zeros(size(xv)), 'x');
            
            s.coord      = V(:,1:2);
            s.connec     = F;
            s.interpType = 'LINEAR';
            mS           = Mesh.create(s);
        end
        
        function computeDirichletDofs(obj)
            coords = obj.globalMesh.coord;
            minX   = min(coords(:,1)); 
            maxX   = max(coords(:,1));
            minY   = min(coords(:,2)); 
            maxY   = max(coords(:,2));
            
            isDirNode = (abs(coords(:,1) - minX) < obj.nodeTol) | ...
                        (abs(coords(:,1) - maxX) < obj.nodeTol) | ...
                        (abs(coords(:,2) - minY) < obj.nodeTol) | ...
                        (abs(coords(:,2) - maxY) < obj.nodeTol);
                        
            obj.dirichletDofs = find(isDirNode);
        end
        
        function computeLHS(obj)
            numSub = prod(obj.numSubdomains);
            kCell  = cell(numSub, 1);
            
            for i = 1:numSub
                localMesh = obj.localMeshes{i};
                uLoc      = LagrangianFunction.create(localMesh, 1, 'P1');
                weakFormK = @(u,v) DP(Grad(v), Grad(u));
                kCell{i}  = IntegrateLHS(weakFormK, uLoc, uLoc, localMesh, 'Domain', 2);
            end
            obj.localStiffness = kCell;
        end
        
        function computeRHS(obj)
            numSub = prod(obj.numSubdomains);
            fCell  = cell(numSub, 1);
            
            for i = 1:numSub
                localMesh = obj.localMeshes{i};
                uLoc      = LagrangianFunction.create(localMesh, 1, 'P1');
                
                weakFormM = @(u,v) DP(v, u);
                mLocal    = IntegrateLHS(weakFormM, uLoc, uLoc, localMesh, 'Domain', 2);
                fLocVec   = mLocal * ones(uLoc.nDofs, 1);
                
                fCell{i}  = fLocVec;
            end
            obj.localForces = fCell;
        end
    end

    % =========================================================
    % PRIVATE METHODS: EXECUTION & SOLVERS
    % =========================================================
    methods (Access = private)
                
        function [uMono, uFeti] = runConvergenceComparison(obj)
            tol = 1e-10;
        
            % CASE 1: Monolithic System Setup and CG Solver 
            uGlobal   = LagrangianFunction.create(obj.globalMesh, 1, 'P1');
            weakFormK = @(u,v) DP(Grad(v), Grad(u));
            kGlobal   = IntegrateLHS(weakFormK, uGlobal, uGlobal, obj.globalMesh, 'Domain', 2);
            
            weakFormM = @(u,v) DP(v, u);
            mGlobal   = IntegrateLHS(weakFormM, uGlobal, uGlobal, obj.globalMesh, 'Domain', 2);
            fGlobal   = mGlobal * ones(uGlobal.nDofs, 1);
            
            freeDofs  = setdiff(1:uGlobal.nDofs, obj.dirichletDofs);
            kFree     = kGlobal(freeDofs, freeDofs);
            fFree     = fGlobal(freeDofs);
            
            uRed  = kFree \ fFree;
            uMono = zeros(uGlobal.nDofs, 1);
            uMono(freeDofs) = uRed;
            
            x0Mono = zeros(length(fFree), 1);
            pId    = @(r) r; 
            
            tic;
            [~, resMono, ~, ~] = PCG.solve(@(x) kFree*x, fFree, x0Mono, pId, tol, uRed);
            timeMono  = toc;
            
            kappaK = 1;
            nDofsMono = length(fFree);
            
            % CASE 2: Unpreconditioned FETI-DP Interface CG 
            dBar   = obj.fetiSolver.assembleProblem();
            x0Feti = zeros(size(dBar));
            lambdaDummy = zeros(size(dBar)); 
            
            tic;
            [~, resFetiDual, ~, ~] = PCG.solve(@(x) obj.fetiSolver.applyFMat(x), dBar, x0Feti, pId, tol, lambdaDummy);
            timeFetiDual = toc;
            
            kappaF = 1; 
            nDofsFeti = length(dBar);
            
            % CASE 3: Preconditioned FETI-DP (Dirichlet) 
            pDir = @(r) obj.fetiSolver.applyDirichletPrecond(r);
            
            tic;
            [lambdaFetiPCG, resFetiDir, ~, ~] = PCG.solve(@(x) obj.fetiSolver.applyFMat(x), dBar, x0Feti, pDir, tol, lambdaDummy);
            timeFetiDir = toc;
            
            kappaPCG = 1;
           
            obj.printComparisonTable(...
                length(resMono), length(resFetiDual), length(resFetiDir), ...
                kappaK, kappaF, kappaPCG, ...
                nDofsMono, nDofsFeti, ...
                prod(obj.numSubdomains), timeMono, timeFetiDual, timeFetiDir);
        
            obj.plotConvergenceComparison(resMono, resFetiDual, resFetiDir, tol);
            
            uFeti = obj.fetiSolver.reconstructGlobalSolution(lambdaFetiPCG, obj.globalMesh.nnodes);
        end
        
    end
    
    % =========================================================
    % PRIVATE METHODS: POST-PROCESSING
    % =========================================================
    methods (Access = private)
        
        function printComparisonTable(~, iter1, iter2, iter3, kK, kF, kPCG, ...
                                      nDofsMono, nDofsFeti, nSub, t1, t2, t3)
            fprintf('\n');
            fprintf('+------------------------------+----------+--------------+----------+-------------+\n');
            fprintf('| Case                         |  Iter.   | kappa (cond) |   DOFs   | CPU Time(s) |\n');
            fprintf('+------------------------------+----------+--------------+----------+-------------+\n');
            fprintf('| Monolithic CG (K*u = f)      | %8d | %12.2f | %8d | %11.4f |\n', iter1, kK,   nDofsMono, t1);
            fprintf('| FETI-DP Interface CG         | %8d | %12.2f | %8d | %11.4f |\n', iter2, kF,   nDofsFeti, t2);
            fprintf('| FETI-DP PCG + Dirichlet      | %8d | %12.2f | %8d | %11.4f |\n', iter3, kPCG, nDofsFeti, t3);
            fprintf('+------------------------------+----------+--------------+----------+-------------+\n');
            fprintf('  Subdomains: %d x %d = %d\n', round(sqrt(nSub)), round(sqrt(nSub)), nSub);
            fprintf('\n');
        end
        
        function plotConvergenceComparison(~, h1, h2, h3, tol)
            figure('Name', 'CG Convergence - Comparison', ...
                   'Color', 'w', 'Position', [100 100 750 480]);
        
            semilogy(1:length(h1), h1, '-o', 'Color', [0.00 0.45 0.74], ...
                     'LineWidth', 1.8, 'MarkerSize', 4, 'DisplayName', 'Monolithic CG');
            hold on;
        
            semilogy(1:length(h2), h2, '-s', 'Color', [0.85 0.33 0.10], ...
                     'LineWidth', 1.8, 'MarkerSize', 4, 'DisplayName', 'FETI-DP Interface CG');
        
            semilogy(1:length(h3), h3, '-^', 'Color', [0.47 0.67 0.19], ...
                     'LineWidth', 1.8, 'MarkerSize', 4, 'DisplayName', 'FETI-DP PCG + Dirichlet Preconditioner');
        
            yline(tol, '--k', 'LineWidth', 1.2, 'Label', sprintf('tol = %.0e', tol), ...
                  'LabelHorizontalAlignment', 'left');
        
            xlabel('Iteration', 'FontSize', 12);
            ylabel('Relative residual', 'FontSize', 12);
            title('CG Convergence - Global Comparison', 'FontSize', 13);
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