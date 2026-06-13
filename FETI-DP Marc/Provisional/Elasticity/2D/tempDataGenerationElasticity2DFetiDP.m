classdef tempDataGenerationElasticity2DFetiDP < handle
    
    properties (Access = public)
        enablePlots       = true
        computeMonolithic = true
        computeKappa      = false
        exportParaview    = false
        outputPrefix      = 'cantilever2d_feti'
        
        numSubdomains     = [10 2]
        nodeTol           = 1e-10
        pcgTol            = 1e-10
        nPerSide          = 5
    end
    
    properties (Access = private)
        globalMesh
        localMeshes
        
        material
        boundaryConditions
        
        localStiffness
        localForces
        
        fetiSolver
    end
    
    methods (Access = public)
        
        % =================================================================
        % CONSTRUCTOR / MAIN
        % =================================================================
        function obj = tempDataGenerationElasticity2DFetiDP()
            if obj.enablePlots
                close all;
            end
            
            % 1. Mesh generation
            referenceMesh = obj.createLatticeMesh();
            s.nsubdomains   = obj.numSubdomains;
            s.meshReference = referenceMesh;
            s.tolSameNode   = obj.nodeTol;
            m = MeshCreatorFromRVE.create(s);
            [obj.globalMesh, obj.localMeshes, ~, ~, ~, ~, ~] = m.create();
            
            % 2. Physics & Boundary Conditions
            obj.material           = obj.createMaterial(obj.globalMesh);
            obj.boundaryConditions = obj.createBoundaryConditions();
            
            % 3. Assembly
            obj.computeLocalMatrices();
            
            % 4. Solver Setup
            obj.fetiSolver = tempFetiDPElasticity( ...
                obj.globalMesh, obj.localMeshes, ...
                obj.localStiffness, obj.localForces, ...
                obj.nodeTol, obj.globalMesh.ndim, ...
                obj.boundaryConditions);
            
            if obj.enablePlots
                obj.fetiSolver.visualizeFetiNodes();
            end
            
            % 5. Solve
            disp('--- FETI-DP (Dirichlet preconditioner) ---');
            tic;
            results = obj.solveFetiDP();
            results.totalWallTime = toc;
            
            % 6. Post-processing
            obj.printResults(results);
            
            if obj.enablePlots
                obj.plotConvergence(results.residual, obj.pcgTol);
                obj.visualizeDeformedMesh(results.uFeti, 1, 'FETI-DP PCG (Cantilever)');
            end
            
            if obj.exportParaview
                if obj.computeMonolithic
                    obj.exportToParaview(results.uMono, 'monolithic');
                end
                obj.exportToParaview(results.uFeti, 'feti_dp');
            end
        end
    end
    
    methods (Access = private)
        
        % =================================================================
        % SOLVER
        % =================================================================
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
                fGlobal = obj.computeGlobalForces(kGlobal);
                results.timeSetupMono = toc;
                
                tic;
                s.stiffness = kGlobal;
                s.forces    = fGlobal;
                [results.uMono, ~] = obj.createProblemSolver().solve(s);
                results.timeSolveMono = toc;
                results.nDofsMono = length(obj.boundaryConditions.free_dofs);
            end
            
            tic;
            [fMat, dBar] = obj.fetiSolver.assembleProblem();
            x0Feti       = zeros(size(dBar));
            lambdaExact  = fMat \ dBar;
            Pdir         = @(r) obj.fetiSolver.applyDirichletPrecond(r);
            results.timeSetupFeti = toc;
            
            tic;
            [lambdaFetiPcg, residual] = PCG.solve(@(x) fMat * x, dBar, x0Feti, Pdir, tol, lambdaExact);
            results.timeSolveFeti = toc;
            
            if obj.computeKappa
                M = obj.fetiSolver.buildPrecondMatrix();
                eigPcg = eig(full(M * fMat));
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
        
        function problemSolver = createProblemSolver(obj)
            s.solverType         = 'REDUCED';
            s.solverMode         = 'DISP';
            s.solver             = DirectSolver();
            s.boundaryConditions = obj.boundaryConditions;
            s.BCApplier          = BCApplier(struct('mesh', obj.globalMesh, 'boundaryConditions', obj.boundaryConditions));
            problemSolver = ProblemSolver(s);
        end
        
        % =================================================================
        % MESH & DOMAIN
        % =================================================================
        
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
        
        function mS = createAuxeticMesh(obj)
            globalLength = 5.0 / obj.numSubdomains(1);
            
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
            % 1. Definir dimensiones objetivo del subdominio
            globalLength = 2.0; 
            globalHeight = 0.4;
            subLength    = globalLength / obj.numSubdomains(1);
            subHeight    = globalHeight / obj.numSubdomains(2);

            % 2. Cargar los datos del archivo .mat de forma robusta
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

            % 5. Creación del objeto Mesh
            s.coord      = coord;
            s.connec     = connec;
            s.interpType = 'LINEAR';
            mS           = Mesh.create(s);
        end
        
        % function mS = createLatticeMesh(obj)
        %     % 1. Cargar los datos del archivo .mat
        %     data = load('mallaLattice.mat');
        % 
        %     % Extraer la variable dinámicamente para evitar errores de nombre
        %     campos = fieldnames(data);        
        %     varName = campos{1};              
        %     meshData = data.(varName);        
        % 
        %     % 2. Extraer coordenadas y conectividades tal cual vienen
        %     % (Ajusta .coord y .connec si tu archivo usa mayúsculas, ej: .COOR y .CN)
        %     coord  = meshData.coord;   
        %     connec = meshData.connec;
        % 
        %     % Si la celosía viene en cuadriláteros (Q4), la dividimos en triángulos (T3)
        %     if size(connec, 2) == 4
        %         connec = [connec(:, [1 2 3]); connec(:, [1 3 4])];
        %     end
        % 
        %     % 3. Creación del objeto Mesh respetando las coordenadas originales
        %     s.coord      = coord;
        %     s.connec     = connec;
        %     s.interpType = 'LINEAR';
        %     mS           = Mesh.create(s);
        % end

        function mult = computeNodeMultiplicity(obj)
            numSub = prod(obj.numSubdomains);
            mult   = zeros(obj.globalMesh.nnodes, 1);
            for i = 1:numSub
                [~, gNodes] = ismembertol(obj.localMeshes{i}.coord, obj.globalMesh.coord, ...
                    obj.nodeTol, 'ByRows', true);
                mult(gNodes) = mult(gNodes) + 1;
            end
        end
        
        % =================================================================
        % PHYSICS & BOUNDARY CONDITIONS
        % =================================================================
        
        function mat = createMaterial(~, mesh)
            eMod   = 200e9;
            nu     = 0.3;
            ePstr  = eMod / (1 - nu^2);
            nuPstr = nu / (1 - nu);
            s.type    = 'ISOTROPIC';
            s.ptype   = 'ELASTIC';
            s.ndim    = mesh.ndim;
            s.young   = ConstantFunction.create(ePstr, mesh);
            s.poisson = ConstantFunction.create(nuPstr, mesh);
            mat       = Material.create(s);
        end
        
        function bc = createBoundaryConditions(obj)
            minX = min(obj.globalMesh.coord(:, 1));
            maxX = max(obj.globalMesh.coord(:, 1));
            isDir = @(coor) abs(coor(:, 1) - minX) < obj.nodeTol;
            isTip = @(coor) abs(coor(:, 1) - maxX) < obj.nodeTol;
            
            sDir.domain    = isDir;
            sDir.direction = [1, 2];
            sDir.value     = 0;
            
            nTipNodes = sum(isTip(obj.globalMesh.coord));
            sPL.domain    = isTip;
            sPL.direction = 2;
            sPL.value     = -500000 / nTipNodes;
            
            s.mesh         = obj.globalMesh;
            s.dirichletFun = DirichletCondition(obj.globalMesh, sDir);
            s.pointloadFun = TractionLoad(obj.globalMesh, sPL, 'DIRAC');
            s.periodicFun  = [];
            bc = BoundaryConditions(s);
        end
        
        % =================================================================
        % MATRIX ASSEMBLY
        % =================================================================
        
        function computeLocalMatrices(obj)
            numSub = prod(obj.numSubdomains);
            kCell  = cell(numSub, 1);
            fCell  = cell(numSub, 1);
            
            nodeMultiplicity = obj.computeNodeMultiplicity();
            fGlobalTraction  = obj.computeGlobalTractionForces();
            dirDofs = obj.boundaryConditions.dirichlet_dofs;
            dirVals = obj.boundaryConditions.dirichlet_vals;
            
            for i = 1:numSub
                localMesh = obj.localMeshes{i};
                uLoc      = LagrangianFunction.create(localMesh, localMesh.ndim, 'P1');
                matLoc    = obj.createMaterial(localMesh);
                weakK     = @(u, v) DDP(SymGrad(v), DDP(matLoc, SymGrad(u)));
                
                kLoc = IntegrateLHS(weakK, uLoc, uLoc, localMesh, 'Domain', 2);
                fLoc = obj.projectGlobalForcesToLocal(fGlobalTraction, localMesh, nodeMultiplicity);
                
                ndim    = localMesh.ndim;
                nnodes  = localMesh.nnodes;
                [~, gNodes] = ismembertol(localMesh.coord, obj.globalMesh.coord, obj.nodeTol, 'ByRows', true);
                
                localToGlobalDofs = reshape(ndim * (gNodes(:)' - 1) + (1:ndim)', [], 1);
                
                [isDir, locIdx] = ismember(localToGlobalDofs, dirDofs);
                dirLocal = find(isDir);
                
                if ~isempty(dirLocal)
                    vals = dirVals(locIdx(isDir));
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
            uGlobal = LagrangianFunction.create(obj.globalMesh, obj.globalMesh.ndim, 'P1');
            weakK   = @(u, v) DDP(SymGrad(v), DDP(obj.material, SymGrad(u)));
            kGlobal = IntegrateLHS(weakK, uGlobal, uGlobal, obj.globalMesh, 'Domain', 2);
        end
        
        function fGlobal = computeGlobalTractionForces(obj)
            uGlobal = LagrangianFunction.create(obj.globalMesh, obj.globalMesh.ndim, 'P1');
            fGlobal = zeros(uGlobal.nDofs, 1);
            tractionLoads = obj.boundaryConditions.tractionFun;
            if isempty(tractionLoads)
                return;
            end
            for k = 1:numel(tractionLoads)
                fGlobal = fGlobal + tractionLoads(k).computeRHS(uGlobal);
            end
        end
        
        function fGlobal = computeGlobalForces(obj, stiffness)
            fGlobal = obj.computeGlobalTractionForces();
            bc      = obj.boundaryConditions;
            dirDofs = bc.dirichlet_dofs;
            dirVals = bc.dirichlet_vals;
            if ~isempty(dirDofs)
                fGlobal = fGlobal - stiffness(:, dirDofs) * dirVals;
            end
        end
        
        function fLoc = projectGlobalForcesToLocal(obj, fGlobal, localMesh, nodeMultiplicity)
            ndim = localMesh.ndim;
            fLoc = zeros(localMesh.nnodes * ndim, 1);
            [~, gNodes] = ismembertol(localMesh.coord, obj.globalMesh.coord, ...
                obj.nodeTol, 'ByRows', true);
            for j = 1:localMesh.nnodes
                mult = nodeMultiplicity(gNodes(j));
                for d = 1:ndim
                    localDof  = (j - 1) * ndim + d;
                    globalDof = (gNodes(j) - 1) * ndim + d;
                    fLoc(localDof) = fGlobal(globalDof) / mult;
                end
            end
        end
        
        % =================================================================
        % POST-PROCESSING & EXPORT 
        % =================================================================
        function printResults(obj, r)
            nSub = prod(obj.numSubdomains);
            
            totalGlobalDofs = obj.globalMesh.nnodes * obj.globalMesh.ndim;
            freeDofs = length(obj.boundaryConditions.free_dofs);
            
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
            title('FETI-DP Preconditioned CG - 2D Elasticity', 'FontSize', 13);
            legend('Location', 'northeast', 'FontSize', 11);
            grid on;
            ax = gca;
            ax.GridAlpha = 0.3;
            ax.YMinorGrid = 'on';
            xlim([1, length(residual) + 1]);
            ylim([tol * 0.1, 2]);
            hold off;
        end
        
        function visualizeDeformedMesh(obj, uGlobal, scaleFactor, titleStr)
            coords   = obj.globalMesh.coord;
            connec   = obj.globalMesh.connec;
            ndim     = obj.globalMesh.ndim;
            numNodes = size(coords, 1);
            uResh     = reshape(uGlobal, ndim, numNodes)';
            defCoords = coords + scaleFactor * uResh;
            dispMag   = sqrt(uResh(:, 1).^2 + uResh(:, 2).^2);
            figure('Name', titleStr, 'Color', 'w');
            hold on; axis equal;
            patch('Faces', connec, 'Vertices', coords, ...
                'FaceColor', 'none', 'EdgeColor', [0.8 0.8 0.8], 'LineStyle', '--');
            patch('Faces', connec, 'Vertices', defCoords, 'FaceVertexCData', dispMag, ...
                'FaceColor', 'interp', 'EdgeColor', '#333333', 'LineWidth', 0.5);
            colormap(jet);
            colorbar;
            title(sprintf('%s (Scale: %gx)', titleStr, scaleFactor));
            xlabel('X'); ylabel('Y');
        end

        function fileBase = exportToParaview(obj, uGlobal, label)
            ndim = obj.globalMesh.ndim;
            numNodes = obj.globalMesh.nnodes;
            
            uResh = reshape(uGlobal, ndim, numNodes)';
            dispMag = sqrt(sum(uResh.^2, 2));
            
            uFun = LagrangianFunction.create(obj.globalMesh, ndim, 'P1');
            uFun.setFValues(uResh);
            
            magFun = LagrangianFunction.create(obj.globalMesh, 1, 'P1');
            magFun.setFValues(dispMag);
            
            fileName = [obj.outputPrefix, '_', label];
            
            s.mesh     = obj.globalMesh;
            s.fun      = {uFun, magFun};
            s.funNames = {'Displacement', 'DisplacementMagnitude'};
            s.type     = 'Paraview';
            s.filename = fileName;
            
            printer = FunctionPrinter.create(s);
            printer.print(); 
            
            fileBase = fullfile(pwd, fileName);
        end       
    end
end