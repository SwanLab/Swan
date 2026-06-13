classdef tempDataGenerationPoisson2DFetiDP < handle
    
    properties (Access = public)
        enablePlots       = true
        computeMonolithic = true
        computeKappa      = true
        numSubdomains     = [2 2]
        nodeTol           = 1e-10
        pcgTol            = 1e-10
        nPerSide          = 3
    end
    
    properties (Access = private)
        globalMesh
        localMeshes
        
        diffusionCoeff
        
        localStiffness
        localForces
        
        fetiSolver
    end
    
    methods (Access = public)
        
        function obj = tempDataGenerationPoisson2DFetiDP()
            if obj.enablePlots
                close all;
            end
            
            referenceMesh = obj.createAlternativeStructuredMesh();
            s.nsubdomains   = obj.numSubdomains;
            s.meshReference = referenceMesh;
            s.tolSameNode   = obj.nodeTol;
            m = MeshCreatorFromRVE.create(s);
            [obj.globalMesh, obj.localMeshes, ~, ~, ~, ~, ~] = m.create();
            
            obj.diffusionCoeff = obj.createDiffusionCoefficient(obj.globalMesh);
            obj.computeLocalMatrices();
            
            obj.fetiSolver = tempFetiDPPoisson( ...
                obj.globalMesh, obj.localMeshes, ...
                obj.localStiffness, obj.localForces, ...
                obj.nodeTol, 1);  % 1 DOF per node for Poisson
            
            if obj.enablePlots
                obj.fetiSolver.visualizeFetiNodes();
            end
            
            disp('--- FETI-DP (Dirichlet preconditioner) ---');
            tic;
            results = obj.solveFetiDP();
            results.totalWallTime = toc;
            
            obj.printResults(results);
            
            if obj.enablePlots
                obj.plotConvergence(results.residual, obj.pcgTol);
                obj.visualizeSolution(results.uFeti, 'FETI-DP PCG (Poisson)');
            end
        end
    end
    
    methods (Access = private)
        
        function results = solveFetiDP(obj)
            tol = obj.pcgTol;
            
            results.computeMonolithic = obj.computeMonolithic;
            results.computeKappa      = obj.computeKappa;
            results.timeSetupMono     = 0;
            results.timeSolveMono     = 0;
            results.nDofsMono         = NaN;
            results.uMono             = [];
            results.relError          = NaN;
            results.kappaPcg          = NaN;
            
            if obj.computeMonolithic
                tic;
                kGlobal = obj.assembleGlobalStiffness();
                fGlobal = obj.computeGlobalForcesWithBC(kGlobal);
                results.timeSetupMono = toc;
                
                coords = obj.globalMesh.coord;
                minX = min(coords(:, 1)); maxX = max(coords(:, 1));
                minY = min(coords(:, 2)); maxY = max(coords(:, 2));
                isOnBoundary = abs(coords(:, 1) - minX) < obj.nodeTol | ...
                              abs(coords(:, 1) - maxX) < obj.nodeTol | ...
                              abs(coords(:, 2) - minY) < obj.nodeTol | ...
                              abs(coords(:, 2) - maxY) < obj.nodeTol;
                dirDofs = find(isOnBoundary);
                freeDofs = setdiff(1:obj.globalMesh.nnodes, dirDofs)';
                
                tic;
                results.uMono = zeros(obj.globalMesh.nnodes, 1);
                kRed = kGlobal(freeDofs, freeDofs);
                fRed = fGlobal(freeDofs);
                results.uMono(freeDofs) = kRed \ fRed;
                results.timeSolveMono = toc;
                results.nDofsMono = length(freeDofs);
            end
            
            tic;
            [fMat, dBar] = obj.fetiSolver.assembleProblem();
            x0Feti       = zeros(size(dBar));
            lambdaExact  = fMat \ dBar;
            pDir         = @(r) obj.fetiSolver.applyDirichletPrecond(r);
            results.timeSetupFeti = toc;
            
            tic;
            [lambdaFetiPcg, residual] = PCG.solve(@(x) fMat * x, dBar, x0Feti, pDir, tol, lambdaExact);
            results.timeSolveFeti = toc;
            
            if obj.computeKappa
                M = obj.fetiSolver.buildPrecondMatrix();
                eigPcg = eig(full(M * fMat));
                results.kappaPcg = max(eigPcg) / min(eigPcg);
            end
            
            results.uFeti     = obj.fetiSolver.reconstructGlobalSolution( ...
                lambdaFetiPcg, obj.globalMesh.nnodes);
            results.nIter     = length(residual);
            results.nDofsFeti = length(dBar);
            results.residual  = residual;
            
            if obj.computeMonolithic
                results.relError = norm(results.uFeti - results.uMono) / norm(results.uMono);
            end
        end
        
        function printResults(obj, r)
            nSub = prod(obj.numSubdomains);
            
            totalGlobalDofs = obj.globalMesh.nnodes;
            coords = obj.globalMesh.coord;
            minX = min(coords(:, 1)); maxX = max(coords(:, 1));
            minY = min(coords(:, 2)); maxY = max(coords(:, 2));
            isOnBoundary = abs(coords(:, 1) - minX) < obj.nodeTol | ...
                          abs(coords(:, 1) - maxX) < obj.nodeTol | ...
                          abs(coords(:, 2) - minY) < obj.nodeTol | ...
                          abs(coords(:, 2) - maxY) < obj.nodeTol;
            freeDofs = sum(~isOnBoundary);
            
            fprintf('\n');
            if r.computeKappa
                fprintf('+---------------------------+-------+--------------+--------+------------------+\n');
                fprintf('| Case                      | Iter. | kappa (cond) |  DOFs  | Total Time (s)   |\n');
                fprintf('+---------------------------+-------+--------------+--------+------------------+\n');
                fprintf('| FETI-DP PCG + Dirichlet   | %5d | %12.2f | %6d | %16.4f |\n', ...
                    r.nIter, r.kappaPcg, r.nDofsFeti, r.timeSetupFeti + r.timeSolveFeti);
                fprintf('+---------------------------+-------+--------------+--------+------------------+\n');
            else
                fprintf('+---------------------------+-------+--------+------------------+\n');
                fprintf('| Case                      | Iter. |  DOFs  | Total Time (s)   |\n');
                fprintf('+---------------------------+-------+--------+------------------+\n');
                fprintf('| FETI-DP PCG + Dirichlet   | %5d | %6d | %16.4f |\n', ...
                    r.nIter, r.nDofsFeti, r.timeSetupFeti + r.timeSolveFeti);
                fprintf('+---------------------------+-------+--------+------------------+\n');
            end
            
            fprintf('  Global DOFs: %d (Total) / %d (Free)\n', totalGlobalDofs, freeDofs);
            fprintf('  Subdomains: %d x %d = %d\n', ...
                round(obj.numSubdomains(1)), round(obj.numSubdomains(2)), nSub);
            fprintf('  Setup (FETI): %.4f s | Solve (PCG): %.4f s\n', ...
                r.timeSetupFeti, r.timeSolveFeti);
                
            if r.computeMonolithic
                fprintf('  Monolithic ref: setup %.4f s | solve %.4f s\n', ...
                    r.timeSetupMono, r.timeSolveMono);
                fprintf('  Relative error (FETI-DP vs monolithic): %e\n', r.relError);
                if r.relError < 1e-10
                    disp('Success: FETI-DP solution matches the monolithic direct solver.');
                end
            end
            fprintf('  Wall-clock (full run): %.4f s\n\n', r.totalWallTime);
        end
        
        function plotConvergence(obj, residual, tol)
            figure('Name', 'FETI-DP PCG Convergence', 'Color', 'w', 'Position', [100 100 650 420]);
            semilogy(1:length(residual), residual, '-^', ...
                'Color', [0.47 0.67 0.19], 'LineWidth', 1.8, 'MarkerSize', 4, ...
                'DisplayName', 'FETI-DP PCG (Dirichlet)');
            hold on;
            yline(tol, '--k', 'LineWidth', 1.2, ...
                'Label', sprintf('tol = %.0e', tol), 'LabelHorizontalAlignment', 'left');
            xlabel('Iteration', 'FontSize', 12);
            ylabel('Relative Residual ||r_k|| / ||r_0||', 'FontSize', 12);
            title('FETI-DP Preconditioned CG - 2D Poisson', 'FontSize', 13);
            legend('Location', 'northeast', 'FontSize', 11);
            grid on;
            ax = gca;
            ax.GridAlpha = 0.3;
            ax.YMinorGrid = 'on';
            xlim([1, length(residual) + 1]);
            ylim([tol * 0.1, 2]);
            hold off;
        end
        
        function computeLocalMatrices(obj)
            numSub = prod(obj.numSubdomains);
            kCell  = cell(numSub, 1);
            fCell  = cell(numSub, 1);
            
            nodeMultiplicity = obj.computeNodeMultiplicity();
            fGlobal          = obj.computeGlobalForceVector();
            
            for i = 1:numSub
                localMesh = obj.localMeshes{i};
                uLoc      = LagrangianFunction.create(localMesh, 1, 'P1');
                weakK = @(u,v) DP(Grad(v), Grad(u));

                
                kLoc = IntegrateLHS(weakK, uLoc, uLoc, localMesh, 'Domain', 2);
                fLoc = obj.projectGlobalForcesToLocal(fGlobal, localMesh, nodeMultiplicity);
                
                [~, gNodes] = ismembertol(localMesh.coord, obj.globalMesh.coord, obj.nodeTol, 'ByRows', true);
                dirDofs = obj.getDirichletDofs(gNodes);
                
                if ~isempty(dirDofs)
                    vals = zeros(length(dirDofs), 1);
                    fLoc = fLoc - kLoc(:, dirDofs) * vals;
                    fLoc(dirDofs) = vals;
                    
                    kLoc(dirDofs, :) = 0;
                    kLoc(:, dirDofs) = 0;
                    kLoc(dirDofs, dirDofs) = speye(length(dirDofs));
                end
                
                kCell{i} = kLoc;
                fCell{i} = fLoc;
            end
            
            obj.localStiffness = kCell;
            obj.localForces    = fCell;
        end
        
        function dirDofs = getDirichletDofs(obj, globalNodes)
            coords = obj.globalMesh.coord(globalNodes, :);
            minX = min(obj.globalMesh.coord(:, 1)); 
            maxX = max(obj.globalMesh.coord(:, 1));
            minY = min(obj.globalMesh.coord(:, 2)); 
            maxY = max(obj.globalMesh.coord(:, 2));
            
            isOnBoundary = abs(coords(:, 1) - minX) < obj.nodeTol | ...
                          abs(coords(:, 1) - maxX) < obj.nodeTol | ...
                          abs(coords(:, 2) - minY) < obj.nodeTol | ...
                          abs(coords(:, 2) - maxY) < obj.nodeTol;
            
            dirDofs = find(isOnBoundary);
        end
        
        function fGlobal = computeGlobalForceVector(obj)
            uGlobal = LagrangianFunction.create(obj.globalMesh, 1, 'P1');
            fGlobal = zeros(uGlobal.nDofs, 1);
            
            sourceFunc = ConstantFunction.create(1, obj.globalMesh);
            weakRHS = @(v) DP(v, sourceFunc);
            fGlobal = IntegrateRHS(weakRHS, uGlobal, obj.globalMesh, 'Domain', 2);
        end
        
        function fLoc = projectGlobalForcesToLocal(obj, fGlobal, localMesh, nodeMultiplicity)
            fLoc = zeros(localMesh.nnodes, 1);
            [~, gNodes] = ismembertol(localMesh.coord, obj.globalMesh.coord, ...
                obj.nodeTol, 'ByRows', true);
            for j = 1:localMesh.nnodes
                mult = nodeMultiplicity(gNodes(j));
                fLoc(j) = fGlobal(gNodes(j)) / mult;
            end
        end
        
        function mult = computeNodeMultiplicity(obj)
            numSub = prod(obj.numSubdomains);
            mult   = zeros(obj.globalMesh.nnodes, 1);
            for i = 1:numSub
                [~, gNodes] = ismembertol(obj.localMeshes{i}.coord, obj.globalMesh.coord, obj.nodeTol, 'ByRows', true);
                mult(gNodes) = mult(gNodes) + 1;
            end
        end
        
        function mS = createStructuredMesh(obj)
            globalLength = 1.0;
            globalHeight = 1.0;
            numSubX      = obj.numSubdomains(1);
            numSubY      = obj.numSubdomains(2);
            subLength    = globalLength / numSubX;
            subHeight    = globalHeight / numSubY;
            x1 = linspace(0, subLength, obj.nPerSide);
            x2 = linspace(0, subHeight, obj.nPerSide);
            [xv, yv] = meshgrid(x1, x2);
            [F, V]   = mesh2tri(xv, yv, zeros(size(xv)), 'x');
            s.coord      = V(:, 1:2);
            s.connec     = F;
            s.interpType = 'LINEAR';
            mS           = Mesh.create(s);
        end

        function mS = createAlternativeStructuredMesh(obj)
            globalLength = 1.0 / obj.numSubdomains(1);
            
            data = load('DEF_Q4auxL_1.mat');
            coord = data.EIFEoper.MESH.COOR;           
            cnQ4 = double(data.EIFEoper.MESH.CN);    
            
            minX = min(coord(:,1));
            maxX = max(coord(:,1));
            minY = min(coord(:,2));
            
            scale = globalLength / (maxX - minX);
            
            coord(:,1) = (coord(:,1) - minX) * scale;
            coord(:,2) = (coord(:,2) - minY) * scale;
            
            s.coord = coord;
            s.connec = [cnQ4(:, [1 2 3]); cnQ4(:, [1 3 4])];
            s.interpType = 'LINEAR';
            mS = Mesh.create(s);
        end
        
        function kFunc = createDiffusionCoefficient(~, mesh)
            kFunc = ConstantFunction.create(1.0, mesh);
        end
        
        function kGlobal = assembleGlobalStiffness(obj)
            uGlobal = LagrangianFunction.create(obj.globalMesh, 1, 'P1');
            weakK = @(u,v) DP(Grad(v), Grad(u));
            kGlobal = IntegrateLHS(weakK, uGlobal, uGlobal, obj.globalMesh, 'Domain', 2);
        end
        
        function fGlobal = computeGlobalForcesWithBC(obj, stiffness)
            fGlobal = obj.computeGlobalForceVector();
            
            coords = obj.globalMesh.coord;
            minX = min(coords(:, 1)); maxX = max(coords(:, 1));
            minY = min(coords(:, 2)); maxY = max(coords(:, 2));
            
            isOnBoundary = abs(coords(:, 1) - minX) < obj.nodeTol | ...
                          abs(coords(:, 1) - maxX) < obj.nodeTol | ...
                          abs(coords(:, 2) - minY) < obj.nodeTol | ...
                          abs(coords(:, 2) - maxY) < obj.nodeTol;
            
            dirDofs = find(isOnBoundary);
            dirVals = zeros(length(dirDofs), 1);
            
            if ~isempty(dirDofs)
                fGlobal = fGlobal - stiffness(:, dirDofs) * dirVals;
            end
        end
        
        function visualizeSolution(obj, uGlobal, titleStr)
            coords   = obj.globalMesh.coord;
            connec   = obj.globalMesh.connec;
            figure('Name', titleStr, 'Color', 'w');
            hold on; axis equal;
            patch('Faces', connec, 'Vertices', coords, 'FaceVertexCData', uGlobal, ...
                'FaceColor', 'interp', 'EdgeColor', '#333333', 'LineWidth', 0.5);
            colormap(jet);
            colorbar;
            title(titleStr);
            xlabel('X'); ylabel('Y');
        end
    end
end
