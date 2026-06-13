classdef tempDataGenerationElasticity3DEdgeAverage < handle
    % Cantilever 3D elasticity benchmark with FETI-DP edge-average primals.
    
    properties (Access = public)
        enablePlots       = false
        computeMonolithic = true
        computeKappa      = false
        
        exportParaview    = true
        outputFolder      = 'paraview_elasticity_3d_edge_average'
        outputPrefix      = 'cantilever3d_edge_average'
        
        numSubdomains     = [10 2 2]
        nElemPerSubdomain = [6 6 6]
        nodeTol           = 1e-10
        pcgTol            = 1e-10
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
        
        function obj = tempDataGenerationElasticity3DEdgeAverage()
            if obj.enablePlots
                close all;
            end
            
            [obj.globalMesh, obj.localMeshes] = obj.createStructuredDecompositionMeshes();
            
            obj.material           = obj.createMaterial(obj.globalMesh);
            obj.boundaryConditions = obj.createBoundaryConditions();
            obj.computeLocalMatrices();
            
            obj.fetiSolver = tempFetiDPElasticity3DEdgeAverage( ...
                obj.globalMesh, obj.localMeshes, ...
                obj.localStiffness, obj.localForces, ...
                obj.nodeTol, obj.globalMesh.ndim, ...
                obj.boundaryConditions);
            
            if obj.enablePlots
                obj.fetiSolver.visualizeFetiNodes();
            end
            
            disp('--- Starting 3D Elasticity FETI-DP (Dirichlet Preconditioner) ---');
            tic;
            results = obj.solveFetiDP();
            results.totalWallTime = toc;
            
            obj.printResults(results);
            
            if obj.exportParaview
                if obj.computeMonolithic
                    monoFile = obj.exportToParaview(results.uMono, 'monolithic');
                    fprintf('  %s.vtu\n', monoFile);
                end
                fetiFile = obj.exportToParaview(results.uFeti, 'feti_dp');
                fprintf('ParaView files written:\n');
                fprintf('  %s.vtu\n', fetiFile);
                fprintf('Open them in ParaView and apply "Warp By Vector" using Displacement.\n\n');
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
            Pdir         = @(r) obj.fetiSolver.applyDirichletPrecond(r);
            results.timeSetupFeti = toc;
            
            lambdaExact = fMat \ dBar;
            
            tic;
            [lambdaFetiPcg, residual] = PCG.solve(@(x) fMat * x, dBar, x0Feti, Pdir, tol, lambdaExact);
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
            fprintf('  Subdomains: %d x %d x %d = %d\n', ...
                round(obj.numSubdomains(1)), round(obj.numSubdomains(2)), round(obj.numSubdomains(3)), nSub);
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
        
        function mS = createStructuredMesh(obj)
            subLength = 5.0 / obj.numSubdomains(1);
            subHeight = 1.0 / obj.numSubdomains(2);
            subWidth  = 1.0 / obj.numSubdomains(3);
            n = obj.nElemPerSubdomain;
            mS = HexaMesh(subLength, subHeight, subWidth, n(1), n(2), n(3));
        end
        
        function [mD, mSbd] = createStructuredDecompositionMeshes(obj)
            referenceMesh = obj.createStructuredMesh();
            refCoord = referenceMesh.coord;
            refConnec = referenceMesh.connec;
            
            refMin = min(refCoord, [], 1);
            refMax = max(refCoord, [], 1);
            subSize = refMax - refMin;
            
            nX = obj.numSubdomains(1);
            nY = obj.numSubdomains(2);
            nZ = obj.numSubdomains(3);
            nSub = prod(obj.numSubdomains);
            
            mSbd = cell(nSub, 1);
            globalCoord = zeros(0, 3);
            globalConnec = zeros(nSub * referenceMesh.nelem, referenceMesh.nnodeElem);
            nextElem = 1;
            
            for k = 1:nZ
                for j = 1:nY
                    for i = 1:nX
                        subId = (k - 1) * nX * nY + (j - 1) * nX + i;
                        offset = [(i - 1) * subSize(1), ...
                                  (j - 1) * subSize(2), ...
                                  (k - 1) * subSize(3)];
                        coord = refCoord + offset;
                        
                        sLoc.coord = coord;
                        sLoc.connec = refConnec;
                        mSbd{subId} = Mesh.create(sLoc);
                        
                        localToGlobal = obj.mapLocalNodesToGlobal(coord, globalCoord);
                        newNodes = localToGlobal == 0;
                        if any(newNodes)
                            nOld = size(globalCoord, 1);
                            nNew = sum(newNodes);
                            globalCoord = [globalCoord; coord(newNodes, :)];
                            localToGlobal(newNodes) = (nOld + 1:nOld + nNew)';
                        end
                        
                        elemRows = nextElem:nextElem + size(refConnec, 1) - 1;
                        globalConnec(elemRows, :) = localToGlobal(refConnec);
                        nextElem = nextElem + size(refConnec, 1);
                    end
                end
            end
            
            sGlob.coord = globalCoord;
            sGlob.connec = globalConnec;
            mD = Mesh.create(sGlob);
        end
        
        function localToGlobal = mapLocalNodesToGlobal(obj, coord, globalCoord)
            nLocal = size(coord, 1);
            localToGlobal = zeros(nLocal, 1);
            if isempty(globalCoord)
                return;
            end
            for inode = 1:nLocal
                diff = abs(globalCoord - coord(inode, :));
                match = find(all(diff < obj.nodeTol, 2), 1);
                if ~isempty(match)
                    localToGlobal(inode) = match;
                end
            end
        end
        
        function bc = createBoundaryConditions(obj)
            minX = min(obj.globalMesh.coord(:, 1));
            maxX = max(obj.globalMesh.coord(:, 1));
            isDir = @(coor) abs(coor(:, 1) - minX) < obj.nodeTol;
            isTip = @(coor) abs(coor(:, 1) - maxX) < obj.nodeTol;
            
            sDir.domain    = isDir;
            sDir.direction = [1, 2, 3];
            sDir.value     = 0;
            
            nTipNodes = sum(isTip(obj.globalMesh.coord));
            sPL.domain    = isTip;
            sPL.direction = 3;
            sPL.value     = -10 / nTipNodes;
            
            s.mesh         = obj.globalMesh;
            s.dirichletFun = DirichletCondition(obj.globalMesh, sDir);
            s.pointloadFun = TractionLoad(obj.globalMesh, sPL, 'DIRAC');
            s.periodicFun  = [];
            bc = BoundaryConditions(s);
        end
        
        function bcApplier = createBCApplier(obj)
            s.mesh = obj.globalMesh;
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
            kCell = cell(numSub, 1);
            fCell = cell(numSub, 1);
            
            nodeMultiplicity = obj.computeNodeMultiplicity();
            fGlobalTraction = obj.computeGlobalTractionForces();

            dirDofs = obj.boundaryConditions.dirichlet_dofs;
            dirVals = obj.boundaryConditions.dirichlet_vals;
            
            for i = 1:numSub
                localMesh = obj.localMeshes{i};
                uLoc = LagrangianFunction.create(localMesh, localMesh.ndim, 'P1');
                matLoc = obj.createMaterial(localMesh);
                weakK = @(u, v) DDP(SymGrad(v), DDP(matLoc, SymGrad(u)));
                
                kLoc = IntegrateLHS(weakK, uLoc, uLoc, localMesh, 'Domain', 2);
                fLoc = obj.projectGlobalForcesToLocal(fGlobalTraction, localMesh, nodeMultiplicity);
                
                ndim = localMesh.ndim;
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
            obj.localForces = fCell;
        end
        
        function mult = computeNodeMultiplicity(obj)
            numSub = prod(obj.numSubdomains);
            mult = zeros(obj.globalMesh.nnodes, 1);
            for i = 1:numSub
                [~, gNodes] = ismembertol(obj.localMeshes{i}.coord, obj.globalMesh.coord, obj.nodeTol, 'ByRows', true);
                mult(gNodes) = mult(gNodes) + 1;
            end
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
                    localDof = (j - 1) * ndim + d;
                    globalDof = (gNodes(j) - 1) * ndim + d;
                    fLoc(localDof) = fGlobal(globalDof) / mult;
                end
            end
        end
        
        function mat = createMaterial(~, mesh)
            eMod = 100000;
            nu = 0.3;
            s.type    = 'ISOTROPIC';
            s.ptype   = 'ELASTIC';
            s.ndim    = mesh.ndim;
            s.young   = ConstantFunction.create(eMod, mesh);
            s.poisson = ConstantFunction.create(nu, mesh);
            mat = Material.create(s);
        end
        
        function kGlobal = assembleGlobalStiffness(obj)
            uGlobal = LagrangianFunction.create(obj.globalMesh, obj.globalMesh.ndim, 'P1');
            weakK = @(u, v) DDP(SymGrad(v), DDP(obj.material, SymGrad(u)));
            kGlobal = IntegrateLHS(weakK, uGlobal, uGlobal, obj.globalMesh, 'Domain', 2);
        end
        
        function fGlobal = computeGlobalForces(obj, stiffness)
            fGlobal = obj.computeGlobalTractionForces();
            dirDofs = obj.boundaryConditions.dirichlet_dofs;
            dirVals = obj.boundaryConditions.dirichlet_vals;
            if ~isempty(dirDofs)
                fGlobal = fGlobal - stiffness(:, dirDofs) * dirVals;
            end
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
            
            s.mesh = obj.globalMesh;
            s.fun = {uFun, magFun};
            s.funNames = {'Displacement', 'DisplacementMagnitude'};
            s.type = 'Paraview';
            s.filename = fileName;
            
            printer = FunctionPrinter.create(s);
            printer.print(); 
            
            if ~exist(obj.outputFolder, 'dir')
                mkdir(obj.outputFolder);
            end
            
            try
                movefile([fileName, '.*'], obj.outputFolder);
            catch
            end
            
            fileBase = fullfile(obj.outputFolder, fileName);
        end
    end
end