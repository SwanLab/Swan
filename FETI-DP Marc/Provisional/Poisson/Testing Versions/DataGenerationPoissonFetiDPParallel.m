classdef DataGenerationPoissonFetiDPParallel < handle    
    
    properties (Access = private)
        % Mesh parameters
        globalMesh       
        localMeshes      
        numSubdomains    
        nodeTol          
        
        % Problem definition
        dirichletDofs   
        localStiffness   
        localForces      
        
        % Solver
        fetiSolver       
    end
    
    methods (Access = public)
        
        function obj = DataGenerationPoissonFetiDPParallel()
            close all;
            
            % 1. Initialization Parameters
            obj.numSubdomains = [25 25]; 
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
            obj.fetiSolver = FetiDPPoissonParallel(obj.globalMesh, obj.localMeshes, ...
                                             obj.localStiffness, obj.localForces, ...
                                             obj.nodeTol, obj.dirichletDofs);
            
            % Visualize domain decomposition nodes
            obj.fetiSolver.visualizeFetiNodes();
            
            % 5. Global Convergence Comparison (3 Cases)
            disp('--- Starting Global Convergence Comparison (3 Cases) ---');
            tic;
            [uMono, uFeti] = obj.runConvergenceComparison();
            totalTime = toc; 
            fprintf('Global comparison time: %.4f seconds\n\n', totalTime);
            
            % 6. Visualization & Error Analysis
            obj.visualizePoissonSolution(uMono, 'Monolithic Solution (Poisson)');
            obj.visualizePoissonSolution(uFeti, 'FETI-DP Solution (Poisson)');
            
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
                nPerSide = 8;
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
            % Identifies nodes on the global outer boundary to apply Dirichlet conditions
            coords = obj.globalMesh.coord;
            minX   = min(coords(:,1)); 
            maxX   = max(coords(:,1));
            minY   = min(coords(:,2)); 
            maxY   = max(coords(:,2));
            
            isDirNode = (abs(coords(:,1) - minX) < obj.nodeTol) | ...
                        (abs(coords(:,1) - maxX) < obj.nodeTol) | ...
                        (abs(coords(:,2) - minY) < obj.nodeTol) | ...
                        (abs(coords(:,2) - maxY) < obj.nodeTol);
                        
            % For Poisson, there is 1 DoF per Node, so node indices match DoF indices
            obj.dirichletDofs = find(isDirNode);
        end
        
        function computeLHS(obj)
            % Assembles local stiffness matrices (K) for each subdomain
            numSub = prod(obj.numSubdomains);
            kCell  = cell(numSub, 1);
            
            parfor i = 1:numSub
                localMesh = obj.localMeshes{i};
                uLoc      = LagrangianFunction.create(localMesh, 1, 'P1');
                weakFormK = @(u,v) DP(Grad(v), Grad(u));
                kCell{i}  = IntegrateLHS(weakFormK, uLoc, uLoc, localMesh, 'Domain', 2);
            end
            obj.localStiffness = kCell;
        end
        
        function computeRHS(obj)
            % Assembles local force vectors (F) for each subdomain
            numSub = prod(obj.numSubdomains);
            fCell  = cell(numSub, 1);
            
            parfor i = 1:numSub
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
            % Runs and compares Monolithic CG, FETI-DP CG, and Preconditioned FETI-DP
            tol = 1e-10;
        
            % CASE 1: Monolithic System Setup and CG Solver 
            uGlobal   = LagrangianFunction.create(obj.globalMesh, 1, 'P1');
            weakFormK = @(u,v) DP(Grad(v), Grad(u));
            kGlobal   = IntegrateLHS(weakFormK, uGlobal, uGlobal, obj.globalMesh, 'Domain', 2);
            
            weakFormM = @(u,v) DP(v, u);
            mGlobal   = IntegrateLHS(weakFormM, uGlobal, uGlobal, obj.globalMesh, 'Domain', 2);
            fGlobal   = mGlobal * ones(uGlobal.nDofs, 1);
            
            % Apply Dirichlet BCs
            freeDofs  = setdiff(1:uGlobal.nDofs, obj.dirichletDofs);
            kFree     = kGlobal(freeDofs, freeDofs);
            fFree     = fGlobal(freeDofs);
            
            % Exact monolithic direct solution (for error tracking)
            uRed  = kFree \ fFree;
            uMono = zeros(uGlobal.nDofs, 1);
            uMono(freeDofs) = uRed;
            
            x0mono = zeros(length(fFree), 1);
            Pid    = @(r) r; % Identity preconditioner (No preconditioner)
            
            tic;
            [~, residual1, ~, ~] = PCG.solve(@(x) kFree*x, fFree, x0mono, Pid, tol, uRed);
            timeMono  = toc;
            kappaK    = condest(kFree);
            nDofsMono = length(fFree);
        
            % CASE 2: Unpreconditioned FETI-DP Interface CG 
            [fMat, dBar] = obj.fetiSolver.assembleProblem();
            x0feti       = zeros(size(dBar));
            
            % Exact multiplier solution (for error tracking in PCG)
            lambdaExact = fMat \ dBar; 
            
            tic;
            [~, residual2, ~, ~] = PCG.solve(@(x) fMat*x, dBar, x0feti, Pid, tol, lambdaExact);
            timeFetiDual = toc;
            kappaF       = condest(fMat);
            nDofsFeti    = length(dBar);
        
            % CASE 3: Preconditioned FETI-DP (Dirichlet) 
            Pdir = @(r) obj.fetiSolver.applyDirichletPrecond(r);
        
            tic;
            [lambdaFetiPCG, residual3, ~, ~] = PCG.solve(@(x) fMat*x, dBar, x0feti, Pdir, tol, lambdaExact);
            timeFetiDir = toc;
            
            % Build dense preconditioner matrix ONLY for condition number analysis
            M = obj.fetiSolver.buildPrecondMatrix();
            kappaPCG = condest(M * fMat);
        
            % Print comparison table and plot residuals
            obj.printComparisonTable(...
                length(residual1), length(residual2), length(residual3), ...
                kappaK, kappaF, kappaPCG, ...
                nDofsMono, nDofsFeti, ...
                prod(obj.numSubdomains), timeMono, timeFetiDual, timeFetiDir);
        
            obj.plotConvergenceComparison(residual1, residual2, residual3, tol);
            
            % Reconstruct global field 'u' using the solved multipliers (lambda)
            uFeti = obj.fetiSolver.reconstructGlobalSolution(lambdaFetiPCG, obj.globalMesh.nnodes);
        end
        
        function uDirect = solveMonolithicPoisson(obj)
            % Standalone utility to solve the global problem using a direct solver
            disp('--- Assembling Global Monolithic System (Poisson) ---');
            uGlobal = LagrangianFunction.create(obj.globalMesh, 1, 'P1');
            
            % Build LHS
            weakFormK = @(u,v) DP(Grad(v), Grad(u));
            kGlobal   = IntegrateLHS(weakFormK, uGlobal, uGlobal, obj.globalMesh, 'Domain', 2);
            
            % Build RHS
            weakFormM = @(u,v) DP(v, u);
            mGlobal   = IntegrateLHS(weakFormM, uGlobal, uGlobal, obj.globalMesh, 'Domain', 2);
            fGlobal   = mGlobal * ones(uGlobal.nDofs, 1);
            
            % Apply Dirichlet DoFs
            freeDofs = setdiff(1:uGlobal.nDofs, obj.dirichletDofs);
            kRed     = kGlobal(freeDofs, freeDofs);
            fRed     = fGlobal(freeDofs);
            
            % Direct inversion
            uRed = kRed \ fRed;
            
            uDirect = zeros(uGlobal.nDofs, 1);
            uDirect(freeDofs) = uRed;
            
            fprintf('Maximum scalar value (Monolithic direct) = %.4e\n', max(abs(uDirect)));
        end
        
    end
    
    % =========================================================
    % PRIVATE METHODS: POST-PROCESSING & VISUALIZATION
    % =========================================================
    methods (Access = private)
        
        function printComparisonTable(~, iter1, iter2, iter3, kK, kF, kPCG, ...
                                      nDofsMono, nDofsFeti, nSub, t1, t2, t3)
            % Prints a strictly formatted ASCII table with the solver results
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
            % Plots the relative residual history for the three CG solves
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

        function visualizePoissonSolution(obj, uGlobal, titleStr)
            % Renders the 2D scalar field mapped onto the global mesh
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
       
    end
end