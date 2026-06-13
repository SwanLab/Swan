classdef tempDataGenerationElasticity2DParaview < handle
    % Cantilever 2D elasticity benchmark with ParaView VTU export.

    properties (Access = public)
        outputFolder = 'paraview_elasticity_2d'
        outputPrefix = 'cantilever'
        deformationScale = 1
    end

    properties (Access = private)
        globalMesh
        localMeshes
        numSubdomains
        nodeTol

        material
        boundaryConditions

        localStiffness
        localForces

        fetiSolver
    end

    methods (Access = public)

        function obj = tempDataGenerationElasticity2DParaview()
            obj.numSubdomains = [5 5];
            obj.nodeTol       = 1e-10;

            referenceMesh = obj.createAlternativeStructuredMesh();
            s.nsubdomains   = obj.numSubdomains;
            s.meshReference = referenceMesh;
            s.tolSameNode   = obj.nodeTol;
            m = MeshCreatorFromRVE.create(s);
            [obj.globalMesh, obj.localMeshes, ~, ~, ~, ~, ~] = m.create();

            obj.material           = obj.createMaterial(obj.globalMesh);
            obj.boundaryConditions = obj.createBoundaryConditions();
            obj.computeLocalMatrices();

            obj.fetiSolver = tempFetiDPElasticity( ...
                obj.globalMesh, obj.localMeshes, ...
                obj.localStiffness, obj.localForces, ...
                obj.nodeTol, obj.globalMesh.ndim, ...
                obj.boundaryConditions);

            disp('--- Starting Global Convergence Comparison (ParaView output) ---');
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
            fprintf('Open them in ParaView and apply "Warp By Vector" using the Displacement field.\n');
        end
    end

    methods (Access = private)

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
            sPL.value     = -10 / nTipNodes;

            dirichletFun = DirichletCondition(obj.globalMesh, sDir);
            pointloadFun = TractionLoad(obj.globalMesh, sPL, 'DIRAC');

            s.mesh         = obj.globalMesh;
            s.dirichletFun = dirichletFun;
            s.pointloadFun = pointloadFun;
            s.periodicFun  = [];
            bc = BoundaryConditions(s);
        end

        function bcApplier = createBCApplier(obj)
            s.mesh               = obj.globalMesh;
            s.boundaryConditions = obj.boundaryConditions;
            bcApplier = BCApplier(s);
        end

        function problemSolver = createProblemSolver(obj)
            s.solverType         = 'REDUCED';
            s.solverMode         = 'DISP';
            s.solver             = DirectSolver();
            s.boundaryConditions = obj.boundaryConditions;
            s.BCApplier          = obj.createBCApplier();
            problemSolver = ProblemSolver(s);
        end

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

                ndim   = localMesh.ndim;
                nnodes = localMesh.nnodes;
                [~, gNodes] = ismembertol(localMesh.coord, obj.globalMesh.coord, obj.nodeTol, 'ByRows', true);

                % localToGlobalDofs = zeros(nnodes * ndim, 1);
                % for j = 1:nnodes
                %     for d = 1:ndim
                %         localToGlobalDofs((j - 1) * ndim + d) = (gNodes(j) - 1) * ndim + d;
                %     end
                % end

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

        function fLoc = projectGlobalForcesToLocal(obj, fGlobal, localMesh, nodeMultiplicity)
            ndim = localMesh.ndim;
            fLoc = zeros(localMesh.nnodes * ndim, 1);
            [~, gNodes] = ismembertol(localMesh.coord, obj.globalMesh.coord, obj.nodeTol, 'ByRows', true);

            for j = 1:localMesh.nnodes
                mult = nodeMultiplicity(gNodes(j));
                for d = 1:ndim
                    localDof  = (j - 1) * ndim + d;
                    globalDof = (gNodes(j) - 1) * ndim + d;
                    fLoc(localDof) = fGlobal(globalDof) / mult;
                end
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
            nPerSide     = 6;
            globalLength = 1.0;
            globalHeight = 1.0;
            numSubX      = obj.numSubdomains(1);
            numSubY      = obj.numSubdomains(2);

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

        function mS = createAlternativeStructuredMesh(obj)
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

        function mat = createMaterial(~, mesh)
            eMod   = 100000;
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

        function kGlobal = assembleGlobalStiffness(obj)
            uGlobal = LagrangianFunction.create(obj.globalMesh, obj.globalMesh.ndim, 'P1');
            weakK   = @(u, v) DDP(SymGrad(v), DDP(obj.material, SymGrad(u)));
            kGlobal = IntegrateLHS(weakK, uGlobal, uGlobal, obj.globalMesh, 'Domain', 2);
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

        function [uMono, uFeti] = runConvergenceComparison(obj)
            tol = 1e-10;

            tic;
            kGlobal  = obj.assembleGlobalStiffness();
            fGlobal  = obj.computeGlobalForces(kGlobal);
            bcApplier = obj.createBCApplier();
            kRed     = bcApplier.fullToReducedMatrixDirichlet(kGlobal);
            fRed     = bcApplier.fullToReducedVectorDirichlet(fGlobal);
            freeDofs = obj.boundaryConditions.free_dofs;
            timeSetupMono = toc;

            tic;
            s.stiffness = kGlobal;
            s.forces    = fGlobal;
            [uMono, ~]  = obj.createProblemSolver().solve(s);
            timeSolveMono = toc;

            x0Mono = zeros(length(fRed), 1);
            pId    = @(r) r;
            uRed   = uMono(freeDofs);
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
            fprintf('| Elasticity Case              | Iter. | kappa (cond) |  DOFs  | Setup Time | Solve Time | Total Time |\n');
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

            ndim     = obj.globalMesh.ndim;
            numNodes = obj.globalMesh.nnodes;
            uResh    = reshape(uGlobal, ndim, numNodes)';
            dispMag  = sqrt(sum(uResh.^2, 2));

            uFun = LagrangianFunction.create(obj.globalMesh, ndim, 'P1');
            uFun.setFValues(uResh);

            magFun = LagrangianFunction.create(obj.globalMesh, 1, 'P1');
            magFun.setFValues(dispMag);

            defFun = LagrangianFunction.create(obj.globalMesh, ndim, 'P1');
            defFun.setFValues(obj.deformationScale * uResh);

            fileBase = fullfile(obj.outputFolder, [obj.outputPrefix, '_', label]);
            s.mesh     = obj.globalMesh;
            s.fun      = {uFun, magFun, defFun};
            s.funNames = {'Displacement', 'DisplacementMagnitude', 'ScaledDisplacement'};
            s.type     = 'Paraview';
            s.filename = fileBase;

            printer = FunctionPrinter.create(s);
            printer.print();
        end
    end
end
