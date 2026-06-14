classdef DataGenerationPoisson3D < handle
    % 3D Poisson benchmark with FETI-DP convergence study.
    
    properties (Access = public)
        useMatrixFree           = true;   
        useEdgeAverage          = true; 
        enablePlots             = true;
        computeMonolithic       = true;
        computeUnpreconditioned = false;
        computeKappa            = false;
        exportParaview          = false;  % Cambiado a true por comodidad de prueba
        outputPrefix            = 'poisson3d_feti';
        
        numSubdomains     = [2 2 2];  % [Nx, Ny, Nz]
        nodeTol           = 1e-8;
        pcgTol            = 1e-8;
        nPerSide          = 5;        
        
        % --- PROPIEDADES PARA LATTICE MESH ---
        useLatticeMesh    = false;                   
        latticeMeshFile   = 'C:\Users\Marc Freixinet\Documents\GitHub\Swan\FETI-DP Marc\Final\mallaLattice3D.mat';   
    end
    
    properties (Access = private)
        globalMesh
        localMeshes
        diffusionCoeff
        boundaryConditions   
        localStiffness
        localForces
        fetiSolver
    end
    
    methods (Access = public)
        % =================================================================
        % CONSTRUCTOR / MAIN
        % =================================================================
        function obj = DataGenerationPoisson3D()
            if obj.enablePlots
                close all;
            end
            
            % 1. Mesh generation
            if obj.useLatticeMesh
                if isfile(obj.latticeMeshFile)
                    loadedData = load(obj.latticeMeshFile);
                    vars = fieldnames(loadedData);
                    
                    if isempty(vars)
                        error('File %s is empty.', obj.latticeMeshFile);
                    end
                    
                    loadedMesh = loadedData.(vars{1});
                    
                    if isstruct(loadedMesh) && isfield(loadedMesh, 'coord') && isfield(loadedMesh, 'connec')
                        sMesh.coord  = loadedMesh.coord;
                        sMesh.connec = loadedMesh.connec;
                        referenceMesh = Mesh.create(sMesh); 
                    else
                        referenceMesh = loadedMesh;
                    end
                    disp('--> Using mallaLattice3D.mat');
                else
                    error('Not Found "%s"', obj.latticeMeshFile);
                end
            else
                referenceMesh = UnitTetraMesh(obj.nPerSide, obj.nPerSide, obj.nPerSide);
                disp('--> Using UnitTetraMesh');
            end
            
            s.nsubdomains   = obj.numSubdomains;
            s.meshReference = referenceMesh;
            s.tolSameNode   = obj.nodeTol;
            m = MeshCreatorFromRVE3D(s);
            [obj.globalMesh, obj.localMeshes, ~, ~, ~, ~, ~] = m.create();
            
            % 2. Physics & Boundary Conditions
            obj.diffusionCoeff = obj.createDiffusionCoefficient(obj.globalMesh);
            obj.boundaryConditions = obj.createBoundaryConditions();
            
            % 3. Assembly
            [nodeMultiplicity, localToGlobalMaps] = obj.computeNodeMultiplicity();
            obj.computeLocalMatrices(nodeMultiplicity, localToGlobalMaps);
            
            % 4. Solver Setup
            dofsPerNode = 1;
            obj.fetiSolver = FetiDPSolver( ...
                obj.globalMesh, obj.localMeshes, ...
                obj.localStiffness, obj.localForces, ...
                obj.nodeTol, dofsPerNode, ...
                obj.boundaryConditions, localToGlobalMaps, ...
                obj.useMatrixFree, obj.useEdgeAverage);
            
            if obj.enablePlots
                obj.fetiSolver.visualizeFetiNodes();
            end
            
            % 5. Solve
            disp('--- FETI-DP 3D ---');
            tic;
            results = obj.solveFetiDP();
            results.totalWallTime = toc;
            
            % 6. Post-processing
            obj.printResults(results);
            if obj.enablePlots
                obj.plotConvergence(results, obj.pcgTol);
                if obj.computeMonolithic
                    obj.visualizeSolution(results.uMono, 'Monolithic Solution (3D Poisson)');
                end
                obj.visualizeSolution(results.uFeti, 'FETI-DP Solution (3D Poisson)');
            end
            
            % --- MODIFICACIÓN: EXPORTACIÓN A PARAVIEW ---
            if obj.exportParaview
                if obj.computeMonolithic
                    obj.exportToParaview(results.uMono, 'monolithic');
                end
                obj.exportToParaview(results.uFeti, 'feti_dp');
                
                % Guarda la malla de referencia (sea Lattice o UnitTetraMesh)
                if obj.useLatticeMesh
                    obj.exportMeshToParaview(referenceMesh, 'reference_lattice');
                else
                    obj.exportMeshToParaview(referenceMesh, 'reference_unit_tetra');
                end
            end
        end
    end
    
    methods (Access = private)
        % =================================================================
        % SOLVER
        % =================================================================
        function results = solveFetiDP(obj)
            tol = obj.pcgTol;
            
            results.computeMonolithic       = obj.computeMonolithic;
            results.computeUnpreconditioned = obj.computeUnpreconditioned;
            results.computeKappa            = obj.computeKappa;
            results.timeSetupMono           = 0;
            results.timeSolveMono           = 0;
            results.timeSolveMonoCG         = NaN;
            results.nDofsMono               = NaN;
            results.uMono                   = [];
            results.relError                = NaN;
            results.kappaK                  = NaN;
            results.kappaF                  = NaN;
            results.kappaPcg                = NaN;
            results.residualMono            = [];
            results.nIterMono               = NaN;
            results.residualUnprec          = [];
            results.nIterUnprec             = NaN;
            results.timeSolveUnprec         = NaN;
            
            if obj.computeMonolithic
                tic;
                kGlobal  = obj.assembleGlobalStiffness();
                fGlobal  = obj.computeGlobalForcesWithBC(kGlobal);
                
                freeDofs = obj.boundaryConditions.free_dofs;
                kRed     = kGlobal(freeDofs, freeDofs);
                fRed     = fGlobal(freeDofs);
                results.timeSetupMono = toc;
                
                tic;
                results.uMono = zeros(obj.globalMesh.nnodes, 1);
                results.uMono(freeDofs) = kRed \ fRed;
                results.timeSolveMono = toc;
                results.nDofsMono     = length(freeDofs);
                
                uRed = results.uMono(freeDofs);
                x0   = zeros(size(fRed));
                Pid  = @(r) r;
                tic;
                [~, residualMono] = PCG.solve(@(x) kRed * x, fRed, x0, Pid, tol, uRed);
                results.timeSolveMonoCG = toc;
                results.residualMono    = residualMono;
                results.nIterMono       = length(residualMono);
                
                if obj.computeKappa
                    eigkRed        = real(eig(full(kRed)));
                    results.kappaK = max(eigkRed) / min(eigkRed);
                end
            end
            
            tic;
            [fMat, dBar] = obj.fetiSolver.assembleProblem();
            x0Feti       = zeros(size(dBar));
            
            if obj.useMatrixFree
                lambdaExact = zeros(size(dBar)); 
            else
                lambdaExact = fMat \ dBar; 
            end
            results.timeSetupFeti = toc;
            
            if isa(fMat, 'function_handle')
                fOperator = fMat;             
            else
                fOperator = @(x) fMat * x;    
            end
            
            if obj.computeUnpreconditioned
                Pid = @(r) r;
                tic;
                [~, residualUnprec] = PCG.solve(fOperator, dBar, x0Feti, Pid, tol, lambdaExact);
                results.timeSolveUnprec = toc;
                results.residualUnprec  = residualUnprec;
                results.nIterUnprec     = length(residualUnprec);
                if obj.computeKappa && ~obj.useMatrixFree
                    eigF           = real(eig(full(fMat)));
                    results.kappaF = max(eigF) / min(eigF);
                end
            end
            
            Pdir = @(r) obj.fetiSolver.applyDirichletPrecond(r);
            tic;
            [lambdaFetiPcg, residual] = PCG.solve(fOperator, dBar, x0Feti, Pdir, tol, lambdaExact);
            results.timeSolveFeti = toc;
            
            if obj.computeKappa
                M = obj.fetiSolver.buildPrecondMatrix();
                if obj.useMatrixFree
                    numDuals = length(dBar); fMatDense = zeros(numDuals);
                    for idx = 1:numDuals
                        eVec = zeros(numDuals, 1); eVec(idx) = 1;
                        fMatDense(:, idx) = fOperator(eVec);
                    end
                    eigPcg = eig(full(M * fMatDense));
                else
                    eigPcg = eig(full(M * fMat));
                end
                results.kappaPcg = max(eigPcg) / min(eigPcg);
            end
            
            results.uFeti     = obj.fetiSolver.reconstructGlobalSolution(lambdaFetiPcg, obj.globalMesh.nnodes);
            results.nIter     = length(residual);
            results.nDofsFeti = length(dBar);
            results.residual  = residual;
            
            if obj.computeMonolithic
                results.relError = norm(results.uFeti - results.uMono) / norm(results.uMono);
            end
        end
        
        % =================================================================
        % NODE MULTIPLICITY
        % =================================================================
        function [mult, localToGlobalMaps] = computeNodeMultiplicity(obj)
            numSub            = prod(obj.numSubdomains);
            mult              = zeros(obj.globalMesh.nnodes, 1);
            localToGlobalMaps = cell(numSub, 1);
            for i = 1:numSub
                [~, gNodes] = ismembertol(obj.localMeshes{i}.coord, obj.globalMesh.coord, ...
                    obj.nodeTol, 'ByRows', true);
                mult(gNodes)         = mult(gNodes) + 1;
                localToGlobalMaps{i} = gNodes;
            end
        end
        
        % =================================================================
        % PHYSICS & BOUNDARY CONDITIONS
        % =================================================================
        function kFunc = createDiffusionCoefficient(~, mesh)
            kFunc = ConstantFunction.create(1.0, mesh);
        end
        
        function f = getBoundaryFunction(obj)
            coords = obj.globalMesh.coord;
            minX = min(coords(:,1)); maxX = max(coords(:,1));
            minY = min(coords(:,2)); maxY = max(coords(:,2));
            minZ = min(coords(:,3)); maxZ = max(coords(:,3));
            tol = obj.nodeTol;
            
            f = @(coor) (abs(coor(:,1) - minX) < tol) | (abs(coor(:,1) - maxX) < tol) | ...
                        (abs(coor(:,2) - minY) < tol) | (abs(coor(:,2) - maxY) < tol) | ...
                        (abs(coor(:,3) - minZ) < tol) | (abs(coor(:,3) - maxZ) < tol);
        end
        
        function bc = createBoundaryConditions(obj)
            isBoundary   = obj.getBoundaryFunction();
            isOnBoundary = isBoundary(obj.globalMesh.coord);
            
            bc.dirichlet_dofs = find(isOnBoundary);
            bc.free_dofs      = setdiff(1:obj.globalMesh.nnodes, bc.dirichlet_dofs)';
            bc.dirichlet_vals = zeros(length(bc.dirichlet_dofs), 1);
        end
        
        % =================================================================
        % MATRIX ASSEMBLY 
        % =================================================================
        function computeLocalMatrices(obj, nodeMultiplicity, localToGlobalMaps)
            numSub = prod(obj.numSubdomains);
            kCell  = cell(numSub, 1);
            fCell  = cell(numSub, 1);
            fGlobal = obj.computeGlobalForceVector();
            
            dirDofsGlobal = obj.boundaryConditions.dirichlet_dofs;
            dirValsGlobal = obj.boundaryConditions.dirichlet_vals;
            
            for i = 1:numSub
                localMesh = obj.localMeshes{i};
                gNodes    = localToGlobalMaps{i};      
                
                uLoc  = LagrangianFunction.create(localMesh, 1, 'P1');
                weakK = @(u, v) DP(Grad(v), Grad(u));
                
                kLoc = IntegrateLHS(weakK, uLoc, uLoc, localMesh, 'Domain', 2);
                fLoc = obj.projectGlobalForcesToLocal(fGlobal, localMesh, nodeMultiplicity, gNodes);
                
                localToGlobalDofs = gNodes(:);
                [isDir, locIdx] = ismember(localToGlobalDofs, dirDofsGlobal);
                dirLocal = find(isDir);
                
                if ~isempty(dirLocal)
                    vals = dirValsGlobal(locIdx(isDir));
                    fLoc = fLoc - kLoc(:, dirLocal) * vals;
                    fLoc(dirLocal) = vals;
                    
                    kLoc(dirLocal, :) = 0;
                    kLoc(:, dirLocal) = 0;
                    kLoc(dirLocal, dirLocal) = speye(length(dirLocal));
                end
                
                kCell{i} = kLoc;
                fCell{i} = fLoc;
            end
            obj.localStiffness = kCell;
            obj.localForces    = fCell;
        end
        
        function kGlobal = assembleGlobalStiffness(obj)
            uGlobal = LagrangianFunction.create(obj.globalMesh, 1, 'P1');
            weakK   = @(u, v) DP(Grad(v), Grad(u));
            kGlobal = IntegrateLHS(weakK, uGlobal, uGlobal, obj.globalMesh, 'Domain', 2);
        end
        
        function fGlobal = computeGlobalForceVector(obj)
            uGlobal    = LagrangianFunction.create(obj.globalMesh, 1, 'P1');
            sourceFunc = ConstantFunction.create(1, obj.globalMesh);
            weakRHS    = @(v) DP(v, sourceFunc);
            fGlobal    = IntegrateRHS(weakRHS, uGlobal, obj.globalMesh, 'Domain', 2);
        end
        
        function fGlobal = computeGlobalForcesWithBC(obj, stiffness)
            fGlobal = obj.computeGlobalForceVector();
            dirDofs = obj.boundaryConditions.dirichlet_dofs;
            dirVals = obj.boundaryConditions.dirichlet_vals;
            
            if ~isempty(dirDofs)
                fGlobal = fGlobal - stiffness(:, dirDofs) * dirVals;
            end
        end
        
        function fLoc = projectGlobalForcesToLocal(~, fGlobal, localMesh, nodeMultiplicity, gNodes)
            fLoc = zeros(localMesh.nnodes, 1);
            for j = 1:localMesh.nnodes
                fLoc(j) = fGlobal(gNodes(j)) / nodeMultiplicity(gNodes(j));
            end
        end
        
        % =================================================================
        % POST-PROCESSING & EXPORT
        % =================================================================
        function printResults(obj, r)
            nSub            = prod(obj.numSubdomains);
            totalGlobalDofs = obj.globalMesh.nnodes;
            freeDofs        = obj.boundaryConditions.free_dofs;
            
            fprintf('\n');
            if r.computeKappa
                fprintf('+------------------------------+-------+--------------+--------+------------+------------+\n');
                fprintf('| Case                         | Iter. | kappa (cond) |  DOFs  | Setup (s)  | Solve (s)  |\n');
                fprintf('+------------------------------+-------+--------------+--------+------------+------------+\n');
                if r.computeMonolithic
                    fprintf('| Monolithic CG (K*u = f)      | %5d | %12.2f | %6d | %10.4f | %10.4f |\n', ...
                        r.nIterMono, r.kappaK, r.nDofsMono, r.timeSetupMono, r.timeSolveMonoCG);
                end
                if r.computeUnpreconditioned
                    fprintf('| FETI-DP dual CG (no prec)    | %5d | %12.2f | %6d | %10.4f | %10.4f |\n', ...
                        r.nIterUnprec, r.kappaF, r.nDofsFeti, r.timeSetupFeti, r.timeSolveUnprec);
                end
                fprintf('| FETI-DP PCG + Dirichlet      | %5d | %12.2f | %6d | %10.4f | %10.4f |\n', ...
                    r.nIter, r.kappaPcg, r.nDofsFeti, r.timeSetupFeti, r.timeSolveFeti);
                fprintf('+------------------------------+-------+--------------+--------+------------+------------+\n');
            else
                fprintf('+------------------------------+-------+--------+------------+------------+\n');
                fprintf('| Case                         | Iter. |  DOFs  | Setup (s)  | Solve (s)  |\n');
                fprintf('+------------------------------+-------+--------+------------+------------+\n');
                if r.computeMonolithic
                    fprintf('| Monolithic CG (K*u = f)      | %5d | %6d | %10.4f | %10.4f |\n', ...
                        r.nIterMono, r.nDofsMono, r.timeSetupMono, r.timeSolveMonoCG);
                end
                if r.computeUnpreconditioned
                    fprintf('| FETI-DP dual CG (no prec)    | %5d | %6d | %10.4f | %10.4f |\n', ...
                        r.nIterUnprec, r.nDofsFeti, r.timeSetupFeti, r.timeSolveUnprec);
                end
                fprintf('| FETI-DP PCG + Dirichlet      | %5d | %6d | %10.4f | %10.4f |\n', ...
                    r.nIter, r.nDofsFeti, r.timeSetupFeti, r.timeSolveFeti);
                fprintf('+------------------------------+-------+--------+------------+------------+\n');
            end
            fprintf('  Global DOFs: %d (Total) / %d (Free)\n', totalGlobalDofs, length(freeDofs));
            fprintf('  Subdomains: %d x %d x %d = %d\n', ...
                obj.numSubdomains(1), obj.numSubdomains(2), obj.numSubdomains(3), nSub);
            
            if r.computeMonolithic
                fprintf('  Monolithic direct solve: %.4f s\n', r.timeSolveMono);
                fprintf('  Relative error (FETI-DP vs monolithic): %e\n', r.relError);
                if r.relError < 1e-8
                    disp('  Success: FETI-DP solution matches the monolithic direct solver.');
                end
            end
            fprintf('  Wall-clock (full run): %.4f s\n\n', r.totalWallTime);
        end
        
        function plotConvergence(~, r, tol)
            figure('Name', 'FETI-DP Convergence Comparison', 'Color', 'w', 'Position', [100 100 750 480]);
            hold on;
            if r.computeMonolithic && ~isempty(r.residualMono)
                semilogy(1:length(r.residualMono), r.residualMono, '-o', ...
                    'Color', [0.00 0.45 0.74], 'LineWidth', 1.8, 'MarkerSize', 4, ...
                    'DisplayName', 'Monolithic CG');
            end
            if r.computeUnpreconditioned && ~isempty(r.residualUnprec)
                semilogy(1:length(r.residualUnprec), r.residualUnprec, '-s', ...
                    'Color', [0.85 0.33 0.10], 'LineWidth', 1.8, 'MarkerSize', 4, ...
                    'DisplayName', 'FETI-DP dual CG (no prec)');
            end
            semilogy(1:length(r.residual), r.residual, '-^', ...
                'Color', [0.47 0.67 0.19], 'LineWidth', 1.8, 'MarkerSize', 4, ...
                'DisplayName', 'FETI-DP PCG (Dirichlet)');
            yline(tol, '--k', 'LineWidth', 1.2, ...
                'Label', sprintf('tol = %.0e', tol), 'LabelHorizontalAlignment', 'left');
            xlabel('Iteration', 'FontSize', 12);
            ylabel('Relative Residual ||r_k|| / ||r_0||', 'FontSize', 12);
            title('CG Convergence - 3D Poisson', 'FontSize', 13);
            legend('Location', 'northeast', 'FontSize', 11);
            grid on;
            ax = gca;
            ax.GridAlpha  = 0.3;
            ax.YMinorGrid = 'on';
            ax.YScale     = 'log';
            allRes = [r.residual(:); r.residualMono(:); r.residualUnprec(:)];
            allRes = allRes(allRes > 0 & isfinite(allRes));
            if ~isempty(allRes)
                ylim([min(allRes) * 0.1, 2]);
            end
            allLengths = [length(r.residual), length(r.residualMono), length(r.residualUnprec)];
            xlim([1, max(allLengths) + 1]);
            hold off;
        end
        
        function visualizeSolution(obj, uGlobal, titleStr)
            coords = obj.globalMesh.coord;
            figure('Name', titleStr, 'Color', 'w');
            hold on; axis equal; view(3); grid on;
            
            scatter3(coords(:,1), coords(:,2), coords(:,3), 40, uGlobal, 'filled');
            colormap(jet); 
            c = colorbar;
            c.Label.String = 'u (solution)';
            title(titleStr);
            xlabel('X'); ylabel('Y'); zlabel('Z');
            hold off;
        end
        
        function exportToParaview(obj, uGlobal, label)
            uFun = LagrangianFunction.create(obj.globalMesh, 1, 'P1');
            uFun.setFValues(uGlobal(:));
            
            fileName = [obj.outputPrefix, '_', label];
            uFun.print(fileName);
        end
        
        % =================================================================
        % NUEVO MÉTODO AUXILIAR PARA EXPORTAR MALLAS DE REFERENCIA
        % =================================================================
        function exportMeshToParaview(obj, mesh, label)
            % Creamos una función Lagrangiana sobre la malla que queremos guardar
            uMesh = LagrangianFunction.create(mesh, 1, 'P1');
            % Inicializamos un vector de ceros del tamaño de nodos de la malla
            uMesh.setFValues(zeros(mesh.nnodes, 1));
            
            fileName = [obj.outputPrefix, '_', label];
            uMesh.print(fileName);
            fprintf('  Malla de referencia exportada a ParaView: %s\n', fileName);
        end
    end
end