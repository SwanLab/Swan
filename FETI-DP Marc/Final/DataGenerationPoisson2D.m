classdef DataGenerationPoisson2D < handle
    % 2D Poisson benchmark with FETI-DP convergence study.
    properties (Access = public)
        useMatrixFree           = true   
        useEdgeAverage          = false
        enablePlots             = true
        computeMonolithic       = true
        computeUnpreconditioned = false
        computeKappa            = false
        exportParaview          = true
        outputPrefix            = 'poisson2d_feti'

        numSubdomains     = [2 2]
        nodeTol           = 1e-12
        pcgTol            = 1e-8
        nPerSide          = 5
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
        function obj = DataGenerationPoisson2D()
            if obj.enablePlots
                close all;
            end
            % 1. Mesh generation
            referenceMesh = obj.createStructuredMesh();

            % --- NUEVO: Exportar la malla de referencia a ParaView ---
            if obj.exportParaview
                % Creamos un campo de desplazamientos (u) ficticio lleno de ceros para la malla de referencia
                uRefFun = LagrangianFunction.create(referenceMesh, referenceMesh.ndim, 'P1');
                uRefFun.setFValues(zeros(referenceMesh.nnodes, referenceMesh.ndim));
                
                % Exportamos el archivo con el prefijo configurado
                fileNameRef = [obj.outputPrefix, '_reference_mesh'];
                uRefFun.print(fileNameRef);
                disp(['Malla de referencia exportada a: ', fileNameRef]);
            end
            % ---------------------------------------------------------

            s.nsubdomains   = obj.numSubdomains;
            s.meshReference = referenceMesh;
            s.tolSameNode   = obj.nodeTol;
            m = MeshCreatorFromRVE.create(s);
            [obj.globalMesh, obj.localMeshes, ~, ~, ~, ~, ~] = m.create();
            
            % 2. Physics
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
            disp('--- FETI-DP ---');
            tic;
            results = obj.solveFetiDP();
            results.totalWallTime = toc;
            
            % 6. Post-processing
            obj.printResults(results);
            if obj.enablePlots
                obj.plotConvergence(results, obj.pcgTol);
                if obj.computeMonolithic
                    obj.visualizeSolution(results.uMono, 'Monolithic Solution (Poisson)');
                end
                obj.visualizeSolution(results.uFeti, 'FETI-DP Solution (Poisson)');
            end

            if obj.exportParaview
                if obj.computeMonolithic
                    obj.exportToParaview(results.uMono, 'monolithic');
                    obj.exportToParaview(results.uFeti, 'feti_dp');
                    obj.exportToParaview(abs(results.uFeti - results.uMono), 'abs_error');
                else
                    obj.exportToParaview(results.uFeti, 'feti_dp');
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
            
            % --- Case 1: Monolithic (direct solve + CG for iteration count)
            if obj.computeMonolithic
                tic;
                kGlobal  = obj.assembleGlobalStiffness();
                fGlobal  = obj.computeGlobalForcesWithBC(kGlobal);
                
                % BCs
                freeDofs = obj.boundaryConditions.free_dofs;
                kRed     = kGlobal(freeDofs, freeDofs);
                fRed     = fGlobal(freeDofs);
                results.timeSetupMono = toc;
                
                % Direct solve (reference solution)
                tic;
                results.uMono = zeros(obj.globalMesh.nnodes, 1);
                results.uMono(freeDofs) = kRed \ fRed;
                results.timeSolveMono = toc;
                results.nDofsMono     = length(freeDofs);
                
                % Unpreconditioned CG on reduced system (for iteration count)
                uRed = results.uMono(freeDofs);
                x0   = zeros(size(fRed));
                Pid  = @(r) r;
                tic;
                [~, residualMono] = PCG.solve(@(x) kRed * x, fRed, x0, Pid, tol, uRed);
                results.timeSolveMonoCG = toc;
                results.residualMono    = residualMono;
                results.nIterMono       = length(residualMono);
                
                if obj.computeKappa
                    % results.kappaK = condest(kRed);
                    eigkRed        = real(eig(full(kRed)));
                    results.kappaK = max(eigkRed) / min(eigkRed);
                end
            end
            
            % --- Assemble FETI-DP interface problem (shared by cases 2 & 3)
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
            
            % --- Case 2: Unpreconditioned FETI-DP dual CG ---------------
            if obj.computeUnpreconditioned
                Pid = @(r) r;
                tic;
                [~, residualUnprec] = PCG.solve(fOperator, dBar, x0Feti, Pid, tol, lambdaExact);
                results.timeSolveUnprec = toc;
                results.residualUnprec  = residualUnprec;
                results.nIterUnprec     = length(residualUnprec);
                if obj.computeKappa
                    eigF           = real(eig(full(fMat)));
                    results.kappaF = max(eigF) / min(eigF);
                end
            end
            
            % --- Case 3: Preconditioned FETI-DP (Dirichlet) -------------
            Pdir = @(r) obj.fetiSolver.applyDirichletPrecond(r);
            tic;
            [lambdaFetiPcg, residual,~,~] = PCG.solve(fOperator, dBar, x0Feti, Pdir, tol, lambdaExact);
            results.timeSolveFeti = toc;

            if obj.computeKappa
                M = obj.fetiSolver.buildPrecondMatrix();
                if obj.useMatrixFree
                    numDuals = length(dBar); fMatDense = zeros(numDuals);
                    for idx = 1:numDuals
                        eVec = zeros(numDuals, 1); eVec(idx) = 1;
                        fMatDense(:, idx) = fMat(eVec);
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
        % MESH & DOMAIN
        % =================================================================
        function mS = createStructuredMesh(obj)
            globalLength = 1;
            globalHeight = 1;
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

        function mS = createAuxeticMesh(obj)
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

        function mS = createLatticeMesh(obj)
            globalLength = 1.0;
            globalHeight = 1.0;
            subLength    = globalLength / obj.numSubdomains(1);
            subHeight    = globalHeight / obj.numSubdomains(2);

            data = load('mallaLattice.mat');
            campos = fieldnames(data);
            varName = campos{1};
            meshData = data.(varName);

            coord  = meshData.coord;
            connec = meshData.connec;

            if size(connec, 2) == 4
                connec = [connec(:, [1 2 3]); connec(:, [1 3 4])];
            end

            minX = min(coord(:,1));
            maxX = max(coord(:,1));
            minY = min(coord(:,2));
            maxY = max(coord(:,2));

            anchoOriginal = maxX - minX;
            altoOriginal  = maxY - minY;

            escalaX = subLength / anchoOriginal;
            escalaY = subHeight / altoOriginal;
            escalaGlobal = min(escalaX, escalaY);

            coord(:,1) = (coord(:,1) - minX) * escalaGlobal;
            coord(:,2) = (coord(:,2) - minY) * escalaGlobal;

            s.coord      = coord;
            s.connec     = connec;
            s.interpType = 'LINEAR';
            mS           = Mesh.create(s);
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
        % PHYSICS
        % =================================================================
        function kFunc = createDiffusionCoefficient(~, mesh)
            kFunc = ConstantFunction.create(1.0, mesh);
        end

        % =================================================================
        % BOUNDARY CONDITIONS
        % =================================================================
        function bc = createBoundaryConditions(obj)
            minX = min(obj.globalMesh.coord(:, 1));
            maxX = max(obj.globalMesh.coord(:, 1));
            minY = min(obj.globalMesh.coord(:, 2));
            maxY = max(obj.globalMesh.coord(:, 2));
            
            isBoundary = @(coor) abs(coor(:, 1) - minX) < obj.nodeTol | ...
                                 abs(coor(:, 1) - maxX) < obj.nodeTol | ...
                                 abs(coor(:, 2) - minY) < obj.nodeTol | ...
                                 abs(coor(:, 2) - maxY) < obj.nodeTol;
            
            sDir.domain    = isBoundary;
            sDir.direction = 1;
            sDir.ndim      = 1; 
            sDir.value     = 0;  
            
            s.mesh         = obj.globalMesh;
            s.ndimf        = 1;  
            s.dirichletFun = DirichletCondition(obj.globalMesh, sDir);
            s.pointloadFun = [];
            s.periodicFun  = [];
            
            bc = BoundaryConditions(s);
        end

        % =================================================================
        % MATRIX ASSEMBLY 
        % =================================================================
        function computeLocalMatrices(obj, nodeMultiplicity, localToGlobalMaps)
            numSub = prod(obj.numSubdomains);
            kCell  = cell(numSub, 1);
            fCell  = cell(numSub, 1);
            fGlobal = obj.computeGlobalForceVector();
            
            dirDofs = obj.boundaryConditions.dirichlet_dofs;
            dirVals = obj.boundaryConditions.dirichlet_vals;
            
            for i = 1:numSub
                localMesh = obj.localMeshes{i};
                gNodes    = localToGlobalMaps{i};      
                uLoc  = LagrangianFunction.create(localMesh, 1, 'P1');
                weakK = @(u, v) DP(Grad(v), Grad(u));
                kLoc = IntegrateLHS(weakK, uLoc, uLoc, localMesh, 'Domain', 2);
                fLoc = obj.projectGlobalForcesToLocal(fGlobal, localMesh, nodeMultiplicity, gNodes);
                
                localToGlobalDofs = gNodes(:); % 1 DoF per node in Poisson
                [isDir, locIdx] = ismember(localToGlobalDofs, dirDofs);
                dirLocal = find(isDir);
                
                if ~isempty(dirLocal)
                    vals = dirVals(locIdx(isDir));   % homogeneous Dirichlet
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
            freeDofs        = obj.boundaryConditions.free_dofs; % Obtenido vía OOP
            
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
            fprintf('  Subdomains: %d x %d = %d\n', ...
                round(obj.numSubdomains(1)), round(obj.numSubdomains(2)), nSub);
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
            title('CG Convergence - 2D Poisson', 'FontSize', 13);
            legend('Location', 'northeast', 'FontSize', 11);
            grid on;
            ax = gca;
            ax.GridAlpha  = 0.3;
            ax.YMinorGrid = 'on';
            ax.YScale     = 'log';
            allRes = [r.residual(:); r.residualMono(:); r.residualUnprec(:)];
            allRes = allRes(allRes > 0 & isfinite(allRes));
            ylim([min(allRes) * 0.1, 2]);
            allLengths = [length(r.residual), length(r.residualMono), length(r.residualUnprec)];
            xlim([1, max(allLengths) + 1]);
            hold off;
        end

        function visualizeSolution(obj, uGlobal, titleStr)
            coords = obj.globalMesh.coord;
            connec = obj.globalMesh.connec;
            figure('Name', titleStr, 'Color', 'w');
            hold on; axis equal;
            patch('Faces', connec, 'Vertices', coords, 'FaceVertexCData', uGlobal, ...
                'FaceColor', 'interp', 'EdgeColor', '#333333', 'LineWidth', 0.5);
            colormap(jet);
            c = colorbar;
            c.Label.String = 'u (solution)';
            title(titleStr);
            xlabel('X'); ylabel('Y');
        end

        function exportToParaview(obj, uGlobal, label)
            uFun = LagrangianFunction.create(obj.globalMesh, 1, 'P1');
            uFun.setFValues(uGlobal(:));
            
            fileName = [obj.outputPrefix, '_', label];
            
            uFun.print(fileName);            
        end
        function Milu = createILUpreconditioner(obj,LHS)
            s.LHS = LHS;
            s.type = 'ILU';
            M = Preconditioner.create(s);
            Milu = @(r) M.apply(r);
        end
    end
    
end