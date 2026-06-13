classdef tempDataGenerationPoisson3DFetiDP < handle
    
    properties (Access = public)
        enablePlots       = true
        computeMonolithic = true
        computeKappa      = true
        numSubdomains     = [2 2 2]   % Número de subdominios [Nx, Ny, Nz]
        nodeTol           = 1e-10
        pcgTol            = 1e-10
        nPerSide          = 5         % Nodos por dimensión para la malla nativa
        outputFolder      = 'paraview_poisson_feti_dp_3d'
        outputPrefix      = 'poisson3d'
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
        
        function obj = tempDataGenerationPoisson3DFetiDP()
            if obj.enablePlots
                close all;
            end
            
            referenceMesh = UnitTetraMesh(obj.nPerSide, obj.nPerSide, obj.nPerSide);
            
            % 2. Generación RVE estructurada
            s.nsubdomains   = obj.numSubdomains;
            s.meshReference = referenceMesh;
            s.tolSameNode   = obj.nodeTol;
            m = MeshCreatorFromRVE3D(s);
            [obj.globalMesh, obj.localMeshes, ~, ~, ~, ~, ~] = m.create();
            
            % 3. Propiedades y Matrices Físicas
            obj.diffusionCoeff = obj.createDiffusionCoefficient(obj.globalMesh);
            obj.computeLocalMatrices();
            
            % 4. Solver FETI-DP (1 DOF por nodo para Poisson)
            obj.fetiSolver = tempFetiDPPoisson3D( ...
                obj.globalMesh, obj.localMeshes, ...
                obj.localStiffness, obj.localForces, ...
                obj.nodeTol, 1);
            
            if obj.enablePlots
                obj.fetiSolver.visualizeFetiNodes();
            end
            
            disp('--- FETI-DP 3D (Dirichlet preconditioner) ---');
            tic;
            results = obj.solveFetiDP();
            results.totalWallTime = toc;
            
            obj.printResults(results);
            
            obj.exportToParaview(results.uFeti, 'feti_dp');
            if obj.computeMonolithic
                obj.exportToParaview(results.uMono, 'monolithic');
            end
            
            if obj.enablePlots
                obj.plotConvergence(results.residual, obj.pcgTol);
                obj.visualizeSolution(results.uFeti, 'FETI-DP PCG Solution (3D Poisson)');
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
            
            % --- RESOLUCIÓN MONOLÍTICA (REFERENCIA) ---
            if obj.computeMonolithic
                tic;
                kGlobal = obj.assembleGlobalStiffness();
                fGlobal = obj.computeGlobalForcesWithBC(kGlobal);
                results.timeSetupMono = toc;
                
                isBoundary = obj.getBoundaryFunction();
                isOnBoundary = isBoundary(obj.globalMesh.coord);
                
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
            
            % --- RESOLUCIÓN FETI-DP ---
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
        
        function computeLocalMatrices(obj)
            numSub = prod(obj.numSubdomains);
            kCell  = cell(numSub, 1);
            fCell  = cell(numSub, 1);
            
            nodeMultiplicity = obj.computeNodeMultiplicity();
            fGlobal          = obj.computeGlobalForceVector();
            
            for i = 1:numSub
                localMesh = obj.localMeshes{i};
                uLoc      = LagrangianFunction.create(localMesh, 1, 'P1');
                weakK     = @(u,v) DP(Grad(v), Grad(u));
                
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
        
        % --- NUEVA LÓGICA DE FRONTERAS ESTILO FUNCIONAL ---
        function f = getBoundaryFunction(obj)
            % Calcula la caja delimitadora del dominio global
            coords = obj.globalMesh.coord;
            minX = min(coords(:,1)); maxX = max(coords(:,1));
            minY = min(coords(:,2)); maxY = max(coords(:,2));
            minZ = min(coords(:,3)); maxZ = max(coords(:,3));
            tol = obj.nodeTol;
            
            % Devuelve una función anónima que evalúa si un nodo está en la frontera
            f = @(coor) (abs(coor(:,1) - minX) < tol) | (abs(coor(:,1) - maxX) < tol) | ...
                        (abs(coor(:,2) - minY) < tol) | (abs(coor(:,2) - maxY) < tol) | ...
                        (abs(coor(:,3) - minZ) < tol) | (abs(coor(:,3) - maxZ) < tol);
        end
        
        function dirDofs = getDirichletDofs(obj, globalNodes)
            coords = obj.globalMesh.coord(globalNodes, :);
            isBoundary = obj.getBoundaryFunction();
            dirDofs = find(isBoundary(coords));
        end
        
        function fGlobal = computeGlobalForceVector(obj)
            uGlobal = LagrangianFunction.create(obj.globalMesh, 1, 'P1');
            sourceFunc = ConstantFunction.create(1, obj.globalMesh);
            weakRHS = @(v) DP(v, sourceFunc);
            fGlobal = IntegrateRHS(weakRHS, uGlobal, obj.globalMesh, 'Domain', 2);
        end
        
        function fLoc = projectGlobalForcesToLocal(obj, fGlobal, localMesh, nodeMultiplicity)
            fLoc = zeros(localMesh.nnodes, 1);
            [~, gNodes] = ismembertol(localMesh.coord, obj.globalMesh.coord, obj.nodeTol, 'ByRows', true);
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
            isBoundary = obj.getBoundaryFunction();
            dirDofs = find(isBoundary(obj.globalMesh.coord));
            
            if ~isempty(dirDofs)
                dirVals = zeros(length(dirDofs), 1);
                fGlobal = fGlobal - stiffness(:, dirDofs) * dirVals;
            end
        end
        
        function printResults(obj, r)
            nSub = prod(obj.numSubdomains);
            totalGlobalDofs = obj.globalMesh.nnodes;
            
            fprintf('\n');
            if r.computeKappa
                fprintf('+---------------------------+-------+--------------+--------+------------------+\n');
                fprintf('| Case                      | Iter. | kappa (cond) |  DOFs  | Total Time (s)   |\n');
                fprintf('+---------------------------+-------+--------------+--------+------------------+\n');
                fprintf('| FETI-DP 3D + Dirichlet    | %5d | %12.2f | %6d | %16.4f |\n', ...
                    r.nIter, r.kappaPcg, r.nDofsFeti, r.timeSetupFeti + r.timeSolveFeti);
                fprintf('+---------------------------+-------+--------------+--------+------------------+\n');
            else
                fprintf('+---------------------------+-------+--------+------------------+\n');
                fprintf('| Case                      | Iter. |  DOFs  | Total Time (s)   |\n');
                fprintf('+---------------------------+-------+--------+------------------+\n');
                fprintf('| FETI-DP 3D + Dirichlet    | %5d | %6d | %16.4f |\n', ...
                    r.nIter, r.nDofsFeti, r.timeSetupFeti + r.timeSolveFeti);
                fprintf('+---------------------------+-------+--------+------------------+\n');
            end
            
            fprintf('  Global DOFs: %d \n', totalGlobalDofs);
            fprintf('  Subdomains: %d x %d x %d = %d\n', ...
                obj.numSubdomains(1), obj.numSubdomains(2), obj.numSubdomains(3), nSub);
            fprintf('  Setup (FETI): %.4f s | Solve (PCG): %.4f s\n', ...
                r.timeSetupFeti, r.timeSolveFeti);
                
            if r.computeMonolithic
                fprintf('  Monolithic ref: setup %.4f s | solve %.4f s\n', ...
                    r.timeSetupMono, r.timeSolveMono);
                fprintf('  Relative error (FETI vs monolithic): %e\n', r.relError);
            end
            fprintf('  Wall-clock (full run): %.4f s\n\n', r.totalWallTime);
        end
        
        function plotConvergence(obj, residual, tol)
            figure('Name', 'FETI-DP PCG Convergence', 'Color', 'w');
            semilogy(1:length(residual), residual, '-^', ...
                'Color', [0.47 0.67 0.19], 'LineWidth', 1.8, 'MarkerSize', 4, ...
                'DisplayName', 'FETI-DP PCG');
            hold on;
            yline(tol, '--k', 'LineWidth', 1.2, ...
                'Label', sprintf('tol = %.0e', tol), 'LabelHorizontalAlignment', 'left');
            xlabel('Iteration', 'FontSize', 12);
            ylabel('Relative Residual ||r_k|| / ||r_0||', 'FontSize', 12);
            title('FETI-DP Preconditioned CG - 3D Poisson', 'FontSize', 13);
            legend('Location', 'northeast'); grid on; hold off;
        end
        
        function visualizeSolution(obj, uGlobal, titleStr)
            coords = obj.globalMesh.coord;
            figure('Name', titleStr, 'Color', 'w');
            hold on; axis equal; view(3); grid on;
            
            scatter3(coords(:,1), coords(:,2), coords(:,3), 40, uGlobal, 'filled');
            colormap(jet); colorbar;
            title(titleStr);
            xlabel('X'); ylabel('Y'); zlabel('Z');
            hold off;
        end

        function fileBase = exportToParaview(obj, uGlobal, label)
            if ~exist(obj.outputFolder, 'dir')
                mkdir(obj.outputFolder);
            end
            uFun = LagrangianFunction.create(obj.globalMesh, 1, 'P1');
            uFun.setFValues(uGlobal(:));
            uFun.print("Poisson3DCube")
            % absFun = LagrangianFunction.create(obj.globalMesh, 1, 'P1');
            % absFun.setFValues(abs(uGlobal(:)));
            % fileBase = fullfile(obj.outputFolder, [obj.outputPrefix, '_', label]);
            % s.mesh = obj.globalMesh;
            % s.fun = {uFun, absFun};
            % s.funNames = {'Solution', 'AbsSolution'};
            % s.type = 'Paraview';
            % s.filename = fileBase;
            % printer = FunctionPrinter.create(s);
            % printer.print();
        end
    end
end