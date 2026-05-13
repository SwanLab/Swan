classdef DataGenerationElasticity2D < handle
    % ELASTICITYFETIDPCONVERGENCESTUDY Main script for data generation, setup,
    % and convergence comparison of the FETI-DP method in 2D Elasticity.
    
    properties (Access = private)
        % Mesh parameters
        globalMesh
        localMeshes
        numSubdomains
        nodeTol
        dofsPerNode
        
        % Problem definition
        dirichletDofs
        localStiffness
        localForces
        
        % Solver
        fetiSolver
    end
    
    % =========================================================
    % PUBLIC METHODS: INITIALIZATION & MAIN EXECUTION
    % =========================================================
    methods (Access = public)
        
        function obj = DataGenerationElasticity2D()
            close all;
            
            % 1. Initialization Parameters
            obj.numSubdomains = [12 8]; 
            obj.nodeTol       = 1e-10;
            obj.dofsPerNode   = 2;
            
            % 2. Global and Local Meshes Setup
            referenceMesh   = obj.createStructuredMesh();
            s.nsubdomains   = obj.numSubdomains;
            s.meshReference = referenceMesh;
            s.tolSameNode   = obj.nodeTol;
            m = MeshCreatorFromRVE.create(s);
            [obj.globalMesh, obj.localMeshes, ~, ~, ~, ~, ~] = m.create();
            
            % 3. Apply Boundary Conditions & Assemble System
            obj.computeDirichletDofs();
            obj.computeLocalMatrices();
            
            % 4. Initialize FETI-DP Solver
            obj.fetiSolver = FetiDPElasticity(obj.globalMesh, obj.localMeshes, ...
                                              obj.localStiffness, obj.localForces, ...
                                              obj.nodeTol, obj.dofsPerNode, obj.dirichletDofs);
            
            % Visualize domain decomposition nodes
            obj.fetiSolver.visualizeFetiNodes();
            
            % 5. Global Convergence Comparison (3 Cases)
            disp('--- Starting Global Convergence Comparison (3 Cases) ---');
            tic;
            [uMono, uFeti] = obj.runConvergenceComparison();
            totalTime = toc; 
            fprintf('Global comparison time: %.4f seconds\n\n', totalTime);
            
            % 6. Visualization & Error Analysis
            obj.visualizeDeformedMesh(uMono, 1, 'Monolithic Solution (Cantilever)');
            obj.visualizeDeformedMesh(uFeti, 1, 'FETI-DP Solution (Cantilever)');
            
            fprintf('Max Y-Displacement (FETI-DP) = %.4e\n', max(abs(uFeti(2:2:end))));
            fprintf('Max Y-Displacement (Direct)  = %.4e\n', max(abs(uMono(2:2:end))));
            
            relError = norm(uFeti - uMono) / norm(uMono);
            fprintf('Relative Error (FETI-DP vs Monolithic): %e\n', relError);
            
            if relError < 1e-10
                disp('Success: FETI-DP solution matches the monolithic direct solver.');
            end
        end
    end
    
    % =========================================================
    % PRIVATE METHODS: SETUP & INITIALIZATION
    % =========================================================
    methods (Access = private)
        
        function computeDirichletDofs(obj)
            % Identifies nodes on the left boundary to apply Dirichlet conditions
            coords = obj.globalMesh.coord;
            minX   = min(coords(:,1));
            
            isDirNode = abs(coords(:,1) - minX) < obj.nodeTol;
            dirNodes  = find(isDirNode);
            
            dDofs = zeros(length(dirNodes) * obj.dofsPerNode, 1);
            for d = 1:obj.dofsPerNode
                dDofs(d:obj.dofsPerNode:end) = (dirNodes - 1) * obj.dofsPerNode + d;
            end
            obj.dirichletDofs = sort(dDofs);
        end
        
        function fGlobal = assembleGlobalForce(obj, totalDofs)
            % Assembles the global force vector with a tip load
            fGlobal   = zeros(totalDofs, 1);
            maxX      = max(obj.globalMesh.coord(:,1));

            tipNodes = find(abs(obj.globalMesh.coord(:,1)-maxX) < obj.nodeTol);
            nTipNodes = length(tipNodes);

            loadValue = -10;
            nodeLoad = loadValue / nTipNodes;
            
            for j = 1:nTipNodes
                yDof = (tipNodes(j) - 1) * obj.dofsPerNode + 2;
                fGlobal(yDof) = fGlobal(yDof) + nodeLoad;                
            end
        end
        
        function computeLocalMatrices(obj)
            % Assembles local stiffness matrices (K) and force vectors (F)
            numSub = prod(obj.numSubdomains);
            kCell  = cell(numSub, 1);
            fCell  = cell(numSub, 1);
            
            maxX      = max(obj.globalMesh.coord(:,1));
            tipNodesGlobal = find(abs(obj.globalMesh.coord(:,1) - maxX) < obj.nodeTol);
            nTipNodes      = length(tipNodesGlobal);
            totalLoad      = -10; 
            nodeLoad       = totalLoad / nTipNodes;
            
            nodeMultiplicity = zeros(obj.globalMesh.nnodes, 1);
            for i = 1:numSub
                [~, gNodes] = ismembertol(obj.localMeshes{i}.coord, obj.globalMesh.coord, obj.nodeTol, 'ByRows', true);
                nodeMultiplicity(gNodes) = nodeMultiplicity(gNodes) + 1;
            end
            
            for i = 1:numSub
                localMesh = obj.localMeshes{i};
                uLoc      = LagrangianFunction.create(localMesh, localMesh.ndim, 'P1');
                matLoc    = obj.createMaterial(localMesh);
                weakK     = @(u,v) DDP(SymGrad(v), DDP(matLoc, SymGrad(u)));
                
                kCell{i}  = IntegrateLHS(weakK, uLoc, uLoc, localMesh, 'Domain', 2);
                fLocVec   = zeros(uLoc.nDofs, 1);
                
                localCoords = localMesh.coord;
                [~, gNodes] = ismembertol(localCoords, obj.globalMesh.coord, obj.nodeTol, 'ByRows', true);
                
                for j = 1:size(localCoords, 1)
                    if abs(localCoords(j,1) - maxX) < obj.nodeTol
                        yDof = (j - 1) * localMesh.ndim + 2;
                        fLocVec(yDof) = fLocVec(yDof) + nodeLoad / nodeMultiplicity(gNodes(j));
                    end
                end
                fCell{i} = fLocVec;
            end
            
            obj.localStiffness = kCell;
            obj.localForces    = fCell;
        end
        
        function mS = createStructuredMesh(obj)
            % Creates the reference mesh for the domain decomposition
            if nargin < 2
                nPerSide = 6; % Número de nodos por lado dentro de CADA subdominio
            end
            globalLength = 1.0; 
            globalHeight = 1.0;
            numSubX = obj.numSubdomains(1);
            numSubY = obj.numSubdomains(2);
            
            subLength = globalLength / numSubX;
            subHeight = globalHeight / numSubY;
            x1 = linspace(0, subLength, nPerSide); 
            x2 = linspace(0, subHeight, nPerSide);
            [xv, yv] = meshgrid(x1, x2);
            [F, V]   = mesh2tri(xv, yv, zeros(size(xv)), 'x');
            
            s.coord      = V(:, 1:2);
            s.connec     = F;
            s.interpType = 'LINEAR';
            mS           = Mesh.create(s);
        end
        
        function mat = createMaterial(obj, mesh)
            % Defines the isotropic elastic material
            [young, poisson] = obj.computeElasticProperties(mesh);
            s.type    = 'ISOTROPIC';
            s.ptype   = 'ELASTIC';
            s.ndim    = mesh.ndim;
            s.young   = young;
            s.poisson = poisson;
            mat       = Material.create(s);
        end
        
        function [young, poisson] = computeElasticProperties(~, mesh)
            % Computes plane stress elastic properties
            eMod   = 100000;
            nu     = 0.3;
            ePstr  = eMod / (1 - nu^2);
            nuPstr = nu   / (1 - nu);
            young   = ConstantFunction.create(ePstr,  mesh);
            poisson = ConstantFunction.create(nuPstr, mesh);
        end
    end

    % =========================================================
    % PRIVATE METHODS: EXECUTION & SOLVERS
    % =========================================================
    methods (Access = private)
        
        function [uMono, uFeti] = runConvergenceComparison(obj)
            % Runs and compares Monolithic CG, FETI-DP CG, and Preconditioned FETI-DP
            tol = 1e-10;
            
            % =========================================================
            % CASE 1: Monolithic System Setup and CG Solver
            % =========================================================
            tic; 
            uGlobal = LagrangianFunction.create(obj.globalMesh, obj.globalMesh.ndim, 'P1');
            mat     = obj.createMaterial(obj.globalMesh);
            weakK   = @(u,v) DDP(SymGrad(v), DDP(mat, SymGrad(u)));
            kGlobal = IntegrateLHS(weakK, uGlobal, uGlobal, obj.globalMesh, 'Domain', 2);
            
            totalDofs = obj.globalMesh.nnodes * obj.dofsPerNode;
            fGlobal   = obj.assembleGlobalForce(totalDofs);
            
            freeDofs = setdiff(1:totalDofs, obj.dirichletDofs);
            kRed     = kGlobal(freeDofs, freeDofs);
            fRed     = fGlobal(freeDofs);
            
            x0mono = zeros(length(fRed), 1);
            Pid    = @(r) r;
            timeSetupMono = toc; 
            
            uRed  = kRed \ fRed;
            uMono = zeros(totalDofs, 1);
            uMono(freeDofs) = uRed;
            
            tic; 
            [~, residual1, ~, ~] = PCG.solve(@(x) kRed * x, fRed, x0mono, Pid, tol, uRed);
            timeSolveMono = toc; 
            
            % lambdaMaxK = eigs(kRed, 1, 'largestabs');
            % lambdaMinK = eigs(kRed, 1, 'smallestabs');
            kappaK     = condest(kRed);%lambdaMaxK / lambdaMinK;
            nDofsMono  = length(fRed);
            
            % =========================================================
            % CASE 2: Unpreconditioned FETI-DP Dual CG
            % =========================================================
            tic; 
            [fMat, dBar] = obj.fetiSolver.assembleProblem();
            x0feti       = zeros(size(dBar));
            timeSetupFetiDual = toc; 
            
            lambdaExact = fMat \ dBar;
            
            tic; 
            [~, residual2, ~, ~] = PCG.solve(@(x) fMat * x, dBar, x0feti, Pid, tol, lambdaExact);
            timeSolveFetiDual = toc; 
            
            E_F       = real(eig(fMat));
            kappaF    = max(E_F) / min(E_F);
            nDofsFeti = length(dBar);
            
            % =========================================================
            % CASE 3: Preconditioned FETI-DP (Dirichlet)
            % =========================================================
            tic; 
            Pdir = @(r) obj.fetiSolver.applyDirichletPrecond(r);
            M    = obj.fetiSolver.buildPrecondMatrix();
            timeSetupFetiDir = toc; 
            
            tic; 
            [lambdaFetiPCG, residual3, ~, ~] = PCG.solve(@(x) fMat * x, dBar, x0feti, Pdir, tol, lambdaExact);
            timeSolveFetiDir = toc; 
            
            E_PCG    = eig(full(M * fMat));
            kappaPCG = max(E_PCG) / min(E_PCG);
            
            % =========================================================
            % Post-Processing
            % =========================================================
            obj.printComparisonTable(...
                length(residual1), length(residual2), length(residual3), ...
                kappaK, kappaF, kappaPCG, ...
                nDofsMono, nDofsFeti, ...
                prod(obj.numSubdomains), ...
                timeSetupMono, timeSolveMono, ...
                timeSetupFetiDual, timeSolveFetiDual, ...
                timeSetupFetiDir, timeSolveFetiDir);
            
            obj.plotConvergenceComparison(residual1, residual2, residual3, tol);
            
            uFeti = obj.fetiSolver.reconstructGlobalSolution(lambdaFetiPCG, obj.globalMesh.nnodes);
        end
    end
    
    % =========================================================
    % PRIVATE METHODS: POST-PROCESSING & VISUALIZATION
    % =========================================================
    methods (Access = private)
        
        function printComparisonTable(obj, iter1, iter2, iter3, kK, kF, kPCG, ...
                                      nDofsMono, nDofsFeti, nSub, ...
                                      tSet1, tSol1, tSet2, tSol2, tSet3, tSol3)
            % Prints a strictly formatted ASCII table with the solver results
            fprintf('\n');
            fprintf('+------------------------------+-------+--------------+--------+------------+------------+------------+\n');
            fprintf('| Elasticity Case              | Iter. | kappa (cond) |  DOFs  | Setup Time | Solve Time | Total Time |\n');
            fprintf('+------------------------------+-------+--------------+--------+------------+------------+------------+\n');
            fprintf('| Monolithic CG (K*u = f)      | %5d | %12.2f | %6d | %8.4f s | %8.4f s | %8.4f s |\n', iter1, kK,   nDofsMono, tSet1, tSol1, tSet1+tSol1);
            fprintf('| FETI-DP Interface CG         | %5d | %12.2f | %6d | %8.4f s | %8.4f s | %8.4f s |\n', iter2, kF,   nDofsFeti, tSet2, tSol2, tSet2+tSol2);
            fprintf('| FETI-DP PCG + Dirichlet      | %5d | %12.2f | %6d | %8.4f s | %8.4f s | %8.4f s |\n', iter3, kPCG, nDofsFeti, tSet2+tSet3, tSol3, tSet2+tSet3+tSol3);
            fprintf('+------------------------------+-------+--------------+--------+------------+------------+------------+\n');
            fprintf('  Subdomains: %d x %d = %d\n', round(obj.numSubdomains(1)), round(obj.numSubdomains(2)), nSub);
            fprintf('\n');
        end
        
        function plotConvergenceComparison(~, h1, h2, h3, tol)
            % Plots the relative residual history for the three CG solves
            figure('Name', 'CG Convergence - Cantilever Beam', 'Color', 'w', 'Position', [100 100 750 480]);
            
            semilogy(1:length(h1), h1, '-o', 'Color', [0.00 0.45 0.74], 'LineWidth', 1.8, 'MarkerSize', 4, 'DisplayName', 'Monolithic CG');
            hold on;
            semilogy(1:length(h2), h2, '-s', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.8, 'MarkerSize', 4, 'DisplayName', 'FETI-DP Dual CG');
            semilogy(1:length(h3), h3, '-^', 'Color', [0.47 0.67 0.19], 'LineWidth', 1.8, 'MarkerSize', 4, 'DisplayName', 'FETI-DP PCG (Dirichlet)');
            
            yline(tol, '--k', 'LineWidth', 1.2, 'Label', sprintf('tol = %.0e', tol), 'LabelHorizontalAlignment', 'left');
            
            xlabel('Iteration', 'FontSize', 12);
            ylabel('Relative Residual ||r_k|| / ||r_0||', 'FontSize', 12);
            title('CG Convergence - 2D Elasticity (Cantilever Beam)', 'FontSize', 13);
            legend('Location', 'northeast', 'FontSize', 11);
            
            grid on;
            ax = gca;              
            ax.GridAlpha = 0.3;     
            ax.YMinorGrid = 'on';   
            xlim([1, max([length(h1), length(h2), length(h3)]) + 1]);
            ylim([tol * 0.1, 2]);
            hold off;
        end

        function visualizeDeformedMesh(obj, uGlobal, scaleFactor, titleStr)
            % Renders the original mesh alongside the deformed mesh
            coords    = obj.globalMesh.coord;
            connec    = obj.globalMesh.connec;
            ndim      = obj.globalMesh.ndim;
            numNodes  = size(coords, 1);
            
            uResh     = reshape(uGlobal, ndim, numNodes)';
            defCoords = coords + scaleFactor * uResh;
            dispMag   = sqrt(uResh(:,1).^2 + uResh(:,2).^2);
            
            figure('Name', titleStr, 'Color', 'w');
            hold on; axis equal;
            
            patch('Faces', connec, 'Vertices', coords, 'FaceColor', 'none', 'EdgeColor', [0.8 0.8 0.8], 'LineStyle', '--');
            patch('Faces', connec, 'Vertices', defCoords, 'FaceVertexCData', dispMag, 'FaceColor', 'interp', 'EdgeColor', '#333333', 'LineWidth', 0.5);
            
            colormap(jet); 
            c = colorbar;
            c.Label.String = '||u||';
            title(sprintf('%s (Scale: %gx)', titleStr, scaleFactor));
            xlabel('X'); ylabel('Y');
        end
    end
end