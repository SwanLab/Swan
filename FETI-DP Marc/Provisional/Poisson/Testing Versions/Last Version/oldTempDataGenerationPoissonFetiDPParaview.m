classdef tempDataGenerationPoissonFetiDPParaview < handle
    % Poisson FETI-DP benchmark with ParaView VTU export.
    % Structured to parallel the generalized elasticity implementation.

    properties (Access = public)
        outputFolder = 'paraview_poisson_feti_dp'
        outputPrefix = 'poisson'
        useAlternativeMesh = false
    end

    properties (Access = private)
        globalMesh
        localMeshes
        numSubdomains
        nodeTol
        pcgTol

        boundaryConditions

        localStiffness
        localForces

        fetiSolver
    end

    methods (Access = public)

        function obj = tempDataGenerationPoissonFetiDPParaview()
            obj.numSubdomains = [5 5];
            obj.nodeTol       = 1e-10;
            obj.pcgTol        = 1e-10;

            if obj.useAlternativeMesh
                referenceMesh = obj.createAlternativeStructuredMesh();
            else
                referenceMesh = obj.createStructuredMesh();
            end

            s.nsubdomains   = obj.numSubdomains;
            s.meshReference = referenceMesh;
            s.tolSameNode   = obj.nodeTol;
            m = MeshCreatorFromRVE.create(s);
            [obj.globalMesh, obj.localMeshes, ~, ~, ~, ~, ~] = m.create();

            obj.boundaryConditions = obj.createBoundaryConditions();
            obj.computeLocalMatrices();

            obj.fetiSolver = tempFetiDPPoissonParaview( ...
                obj.globalMesh, obj.localMeshes, ...
                obj.localStiffness, obj.localForces, ...
                obj.nodeTol, 1, ... % dofsPerNode for Poisson
                obj.boundaryConditions);

            disp('--- Starting Poisson FETI-DP comparison (ParaView output) ---');
            tic;
            [uMono, uFeti] = obj.runConvergenceComparison();
            fprintf('Global comparison time: %.4f seconds\n\n', toc);

            monoFile = obj.exportToParaview(uMono, 'monolithic');
            fetiFile = obj.exportToParaview(uFeti, 'feti_dp');

            relError = norm(uFeti - uMono) / norm(uMono);
            fprintf('Relative Error (FETI-DP vs Monolithic): %e\n', relError);
            if relError < 1e-10
                disp('Success: FETI-DP solution matches the monolithic direct solver.');
            end

            fprintf('ParaView files written:\n');
            fprintf('  %s.vtu\n', monoFile);
            fprintf('  %s.vtu\n', fetiFile);
            fprintf('Open them in ParaView and color by Solution or AbsSolution.\n');
        end
    end

    methods (Access = private)

        function bc = createBoundaryConditions(obj)
            minX = min(obj.globalMesh.coord(:, 1));
            maxX = max(obj.globalMesh.coord(:, 1));
            minY = min(obj.globalMesh.coord(:, 2));
            maxY = max(obj.globalMesh.coord(:, 2));

            isDir = @(coor) (abs(coor(:, 1) - minX) < obj.nodeTol) | ...
                            (abs(coor(:, 1) - maxX) < obj.nodeTol) | ...
                            (abs(coor(:, 2) - minY) < obj.nodeTol) | ...
                            (abs(coor(:, 2) - maxY) < obj.nodeTol);

            sDir.domain    = isDir;
            sDir.direction = 1;
            sDir.value     = 0;
            sDir.ndim      = 1;

            dirichletFun = DirichletCondition(obj.globalMesh, sDir);

            s.mesh         = obj.globalMesh;
            s.dirichletFun = dirichletFun;
            s.pointloadFun = [];
            s.periodicFun  = [];
            bc = BoundaryConditions(s);
        end

        function bcApplier = createBCApplier(obj)
            s.mesh               = obj.globalMesh;
            s.boundaryConditions = obj.boundaryConditions;
            bcApplier = BCApplier(s);
        end

        function computeLocalMatrices(obj)
            numSub = prod(obj.numSubdomains);
            kCell  = cell(numSub, 1);
            fCell  = cell(numSub, 1);

            for i = 1:numSub
                localMesh = obj.localMeshes{i};
                uLoc      = LagrangianFunction.create(localMesh, 1, 'P1');

                weakFormK = @(u, v) DP(Grad(v), Grad(u));
                weakFormM = @(u, v) DP(v, u);

                kLoc = IntegrateLHS(weakFormK, uLoc, uLoc, localMesh, 'Domain', 2);
                mLoc = IntegrateLHS(weakFormM, uLoc, uLoc, localMesh, 'Domain', 2);
                fLoc = mLoc * ones(uLoc.nDofs, 1);

                dirDofs = obj.boundaryConditions.dirichlet_dofs;
                dirVals = obj.boundaryConditions.dirichlet_vals;

                [~, gNodes] = ismembertol(localMesh.coord, obj.globalMesh.coord, obj.nodeTol, 'ByRows', true);
                
                % For scalar Poisson, local DOFs directly correspond to node indices
                [isDir, locIdx] = ismember(gNodes, dirDofs);
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
            uGlobal   = LagrangianFunction.create(obj.globalMesh, 1, 'P1');
            weakFormK = @(u, v) DP(Grad(v), Grad(u));
            kGlobal   = IntegrateLHS(weakFormK, uGlobal, uGlobal, obj.globalMesh, 'Domain', 2);
        end

        function fGlobal = computeGlobalForces(obj, stiffness)
            uGlobal   = LagrangianFunction.create(obj.globalMesh, 1, 'P1');
            weakFormM = @(u, v) DP(v, u);
            mGlobal   = IntegrateLHS(weakFormM, uGlobal, uGlobal, obj.globalMesh, 'Domain', 2);
            fGlobal   = mGlobal * ones(uGlobal.nDofs, 1);

            bc      = obj.boundaryConditions;
            dirDofs = bc.dirichlet_dofs;
            dirVals = bc.dirichlet_vals;
            if ~isempty(dirDofs)
                fGlobal = fGlobal - stiffness(:, dirDofs) * dirVals;
            end
        end

        function [uMono, uFeti] = runConvergenceComparison(obj)
            tol = obj.pcgTol;

            tic;
            kGlobal   = obj.assembleGlobalStiffness();
            fGlobal   = obj.computeGlobalForces(kGlobal);
            bcApplier = obj.createBCApplier();
            kRed      = bcApplier.fullToReducedMatrixDirichlet(kGlobal);
            fRed      = bcApplier.fullToReducedVectorDirichlet(fGlobal);
            freeDofs  = obj.boundaryConditions.free_dofs;
            dirDofs   = obj.boundaryConditions.dirichlet_dofs;
            dirVals   = obj.boundaryConditions.dirichlet_vals;
            timeSetupMono = toc;

            tic;
            uRed = kRed \ fRed;
            uMono = zeros(obj.globalMesh.nnodes, 1);
            uMono(freeDofs) = uRed;
            uMono(dirDofs)  = dirVals;
            timeSolveMono = toc;

            x0Mono = zeros(length(fRed), 1);
            pId    = @(r) r;
            [~, residual1, ~, ~] = PCG.solve(@(x) kRed * x, fRed, x0Mono, pId, tol, uRed);
            kappaK    = condest(kRed);
            nDofsMono = length(fRed);

            tic;
            [fMat, dBar] = obj.fetiSolver.assembleProblem();
            x0Feti       = zeros(size(dBar));
            timeSetupFetiDual = toc;

            lambdaExact = fMat \ dBar;

            tic;
            [~, residual2, ~, ~] = PCG.solve(@(x) fMat * x, dBar, x0Feti, pId, tol, lambdaExact);
            timeSolveFetiDual = toc;

            eigF      = real(eig(full(fMat)));
            kappaF    = max(eigF) / min(eigF);
            nDofsFeti = length(dBar);

            tic;
            pDir = @(r) obj.fetiSolver.applyDirichletPrecond(r);
            M    = obj.fetiSolver.buildPrecondMatrix();
            timeSetupFetiDir = toc;

            tic;
            [lambdaFetiPcg, residual3, ~, ~] = PCG.solve(@(x) fMat * x, dBar, x0Feti, pDir, tol, lambdaExact);
            timeSolveFetiDir = toc;

            eigPcg   = eig(full(M * fMat));
            kappaPcg = max(eigPcg) / min(eigPcg);

            obj.printComparisonTable( ...
                length(residual1), length(residual2), length(residual3), ...
                kappaK, kappaF, kappaPcg, ...
                nDofsMono, nDofsFeti, prod(obj.numSubdomains), ...
                timeSetupMono, timeSolveMono, ...
                timeSetupFetiDual, timeSolveFetiDual, ...
                timeSetupFetiDir, timeSolveFetiDir);

            uFeti = obj.fetiSolver.reconstructGlobalSolution(lambdaFetiPcg, obj.globalMesh.nnodes);
        end

        function printComparisonTable(obj, iter1, iter2, iter3, kK, kF, kPCG, ...
                nDofsMono, nDofsFeti, nSub, tSet1, tSol1, tSet2, tSol2, tSet3, tSol3)
            fprintf('\n');
            fprintf('+------------------------------+-------+--------------+--------+------------+------------+------------+\n');
            fprintf('| Poisson Case                 | Iter. | kappa (cond) |  DOFs  | Setup Time | Solve Time | Total Time |\n');
            fprintf('+------------------------------+-------+--------------+--------+------------+------------+------------+\n');
            fprintf('| Monolithic CG (K*u = f)      | %5d | %12.2f | %6d | %8.4f s | %8.4f s | %8.4f s |\n', iter1, kK, nDofsMono, tSet1, tSol1, tSet1 + tSol1);
            fprintf('| FETI-DP Interface CG         | %5d | %12.2f | %6d | %8.4f s | %8.4f s | %8.4f s |\n', iter2, kF, nDofsFeti, tSet2, tSol2, tSet2 + tSol2);
            fprintf('| FETI-DP PCG + Dirichlet      | %5d | %12.2f | %6d | %8.4f s | %8.4f s | %8.4f s |\n', iter3, kPCG, nDofsFeti, tSet2 + tSet3, tSol3, tSet2 + tSet3 + tSol3);
            fprintf('+------------------------------+-------+--------------+--------+------------+------------+------------+\n');
            fprintf('  Subdomains: %d x %d = %d\n', round(obj.numSubdomains(1)), round(obj.numSubdomains(2)), nSub);
            fprintf('\n');
        end

        function fileBase = exportToParaview(obj, uGlobal, label)
            if ~exist(obj.outputFolder, 'dir')
                mkdir(obj.outputFolder);
            end

            uFun = LagrangianFunction.create(obj.globalMesh, 1, 'P1');
            uFun.setFValues(uGlobal(:));

            absFun = LagrangianFunction.create(obj.globalMesh, 1, 'P1');
            absFun.setFValues(abs(uGlobal(:)));

            fileBase = fullfile(obj.outputFolder, [obj.outputPrefix, '_', label]);
            s.mesh = obj.globalMesh;
            s.fun = {uFun, absFun};
            s.funNames = {'Solution', 'AbsSolution'};
            s.type = 'Paraview';
            s.filename = fileBase;

            printer = FunctionPrinter.create(s);
            printer.print();
        end

        function mS = createStructuredMesh(obj)
            nPerSide = 5;
            x1 = linspace(0, 1 / obj.numSubdomains(1), nPerSide);
            x2 = linspace(0, 1 / obj.numSubdomains(2), nPerSide);
            [xv, yv] = meshgrid(x1, x2);
            [F, V] = mesh2tri(xv, yv, zeros(size(xv)), 'x');

            s.coord = V(:, 1:2);
            s.connec = F;
            s.interpType = 'LINEAR';
            mS = Mesh.create(s);
        end

        function mS = createAlternativeStructuredMesh(obj)
            globalLength = 5.0 / obj.numSubdomains(1);

            data = load('DEF_Q4auxL_1.mat');
            coord = data.EIFEoper.MESH.COOR;
            cnQ4 = double(data.EIFEoper.MESH.CN);

            minX = min(coord(:, 1));
            maxX = max(coord(:, 1));
            minY = min(coord(:, 2));

            scale = globalLength / (maxX - minX);
            coord(:, 1) = (coord(:, 1) - minX) * scale;
            coord(:, 2) = (coord(:, 2) - minY) * scale;

            s.coord = coord;
            s.connec = [cnQ4(:, [1 2 3]); cnQ4(:, [1 3 4])];
            s.interpType = 'LINEAR';
            mS = Mesh.create(s);
        end
    end
end
