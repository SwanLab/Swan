classdef DataGenerationElasticityNofMat < handle
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
        function obj = DataGenerationElasticityNofMat()
            close all;
            
            obj.numSubdomains = [70 70];
            obj.nodeTol = 1e-10;
            
            refMesh = obj.createStructuredMesh();
            s.nsubdomains = obj.numSubdomains;
            s.meshReference = refMesh;
            s.tolSameNode = obj.nodeTol;
            
            m = MeshCreatorFromRVE.create(s);
            [obj.globalMesh, obj.localMeshes, ~, ~, ~, ~, ~] = m.create();
            
            obj.material = obj.createMaterial(obj.globalMesh);
            obj.boundaryConditions = obj.createBoundaryConditions();
            
            obj.computeLocalMatrices();
            
            obj.fetiSolver = FetiDPElasticityNofMat( ...
                obj.globalMesh, obj.localMeshes, ...
                obj.localStiffness, obj.localForces, ...
                obj.nodeTol, obj.globalMesh.ndim);
            
            disp('--- Starting Matrix-Free FETI-DP PCG Solver ---');
            uFeti = obj.runFetiSolver();
            
            obj.visualizeDeformedMesh(uFeti, 1, 'FETI-DP Solution');
        end
    end
    
    methods (Access = private)
        function uFeti = runFetiSolver(obj)
            tol = 1e-10;
            
            tic;
            dBar = obj.fetiSolver.setupProblem();
            timeSetup = toc;
            
            x0Feti = zeros(size(dBar));
            dummyExact = zeros(size(dBar));
            
            aFun = @(x) obj.fetiSolver.applyFMat(x);
            pDir = @(r) obj.fetiSolver.applyPreconditioner(r);
            
            tic;
            [lambdaFeti, res, ~, ~] = PCG.solve(aFun, dBar, x0Feti, pDir, tol, dummyExact);
            timeSolve = toc;
            
            uFeti = obj.fetiSolver.reconstructGlobalSolution(lambdaFeti, obj.globalMesh.nnodes);
            
            globalDoFs = obj.globalMesh.nnodes * obj.globalMesh.ndim;
            
            fprintf('Global DoFs:           %d\n', globalDoFs);
            fprintf('Dual (Interface) DoFs: %d\n', length(dBar));
            fprintf('FETI-DP Setup Time:    %.4f s\n', timeSetup);
            fprintf('FETI-DP Solve Time:    %.4f s\n', timeSolve);
            fprintf('Total Time:            %.4f s\n', timeSetup + timeSolve);
            fprintf('PCG Iterations:        %d\n\n', length(res));
            
            obj.plotConvergence(res, tol);
        end
        
        function bc = createBoundaryConditions(obj)
            minX = min(obj.globalMesh.coord(:,1));
            maxX = max(obj.globalMesh.coord(:,1));
            
            sDir.domain = @(coor) abs(coor(:,1) - minX) < obj.nodeTol;
            sDir.direction = [1, 2];
            sDir.value = 0;
            
            isTip = @(coor) abs(coor(:,1) - maxX) < obj.nodeTol;
            nTipNodes = sum(isTip(obj.globalMesh.coord));
            
            sPL.domain = isTip;
            sPL.direction = 2;
            sPL.value = -10 / nTipNodes;
            
            s.mesh = obj.globalMesh;
            s.dirichletFun = DirichletCondition(obj.globalMesh, sDir);
            s.pointloadFun = TractionLoad(obj.globalMesh, sPL, 'DIRAC');
            s.periodicFun = [];
            
            bc = BoundaryConditions(s);
        end
        
        function computeLocalMatrices(obj)
            numSub = prod(obj.numSubdomains);
            kCell = cell(numSub, 1);
            fCell = cell(numSub, 1);
            
            mult = obj.computeNodeMultiplicity();
            fTraction = obj.computeGlobalTractionForces();
            
            dirDofs = obj.boundaryConditions.dirichlet_dofs;
            dirVals = obj.boundaryConditions.dirichlet_vals;
            
            for i = 1:numSub
                localMesh = obj.localMeshes{i};
                uLoc = LagrangianFunction.create(localMesh, localMesh.ndim, 'P1');
                matLoc = obj.createMaterial(localMesh);
                weakK = @(u, v) DDP(SymGrad(v), DDP(matLoc, SymGrad(u)));
                
                kLoc = IntegrateLHS(weakK, uLoc, uLoc, localMesh, 'Domain', 2);
                fLoc = obj.projectGlobalForcesToLocal(fTraction, localMesh, mult);
           
                ndim = localMesh.ndim;
                nnodes = localMesh.nnodes;
                [~, gNodes] = ismembertol(localMesh.coord, obj.globalMesh.coord, obj.nodeTol, 'ByRows', true);
                
                locToGlob = zeros(nnodes * ndim, 1);
                for j = 1:nnodes
                    for d = 1:ndim
                        locToGlob((j - 1) * ndim + d) = (gNodes(j) - 1) * ndim + d;
                    end
                end
                
                [isDir, locIdx] = ismember(locToGlob, dirDofs);
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
        
        function fGlobal = computeGlobalTractionForces(obj)
            uGlobal = LagrangianFunction.create(obj.globalMesh, obj.globalMesh.ndim, 'P1');
            fGlobal = zeros(uGlobal.nDofs, 1);
            loads = obj.boundaryConditions.tractionFun;
            if isempty(loads); return; end
            for k = 1:numel(loads)
                fGlobal = fGlobal + loads(k).computeRHS(uGlobal);
            end
        end
        
        function fLoc = projectGlobalForcesToLocal(obj, fGlobal, localMesh, mult)
            ndim = localMesh.ndim;
            fLoc = zeros(localMesh.nnodes * ndim, 1);
            [~, gNodes] = ismembertol(localMesh.coord, obj.globalMesh.coord, obj.nodeTol, 'ByRows', true);
            for j = 1:localMesh.nnodes
                for d = 1:ndim
                    fLoc((j - 1) * ndim + d) = fGlobal((gNodes(j) - 1) * ndim + d) / mult(gNodes(j));
                end
            end
        end
        
        function mult = computeNodeMultiplicity(obj)
            numSub = prod(obj.numSubdomains);
            mult = zeros(obj.globalMesh.nnodes, 1);
            for i = 1:numSub
                [~, gNodes] = ismembertol(obj.localMeshes{i}.coord, obj.globalMesh.coord, obj.nodeTol, 'ByRows', true);
                mult(gNodes) = mult(gNodes) + 1;
            end
        end
        
        function mS = createStructuredMesh(obj)
            nPerSide = 6;
            x1 = linspace(0, 1.0 / obj.numSubdomains(1), nPerSide);
            x2 = linspace(0, 1.0 / obj.numSubdomains(2), nPerSide);
            [xv, yv] = meshgrid(x1, x2);
            [F, V] = mesh2tri(xv, yv, zeros(size(xv)), 'x');
            
            s.coord = V(:, 1:2);
            s.connec = F;
            s.interpType = 'LINEAR';
            mS = Mesh.create(s);
        end
        
        function mat = createMaterial(~, mesh)
            eMod = 100000;
            nu = 0.3;
            s.type = 'ISOTROPIC';
            s.ptype = 'ELASTIC';
            s.ndim = mesh.ndim;
            s.young = ConstantFunction.create(eMod / (1 - nu^2), mesh);
            s.poisson = ConstantFunction.create(nu / (1 - nu), mesh);
            mat = Material.create(s);
        end
        
        function plotConvergence(~, res, tol)
            figure('Name', 'PCG Convergence', 'Color', 'w');
            semilogy(1:length(res), res, '-^', 'Color', [0.47 0.67 0.19], 'LineWidth', 1.8);
            yline(tol, '--k', 'LineWidth', 1.2, 'Label', sprintf('tol = %.0e', tol));
            xlabel('Iteration'); ylabel('Relative Residual');
            title('PCG Convergence (Matrix-Free FETI-DP)');
            grid on;
        end
        
        function visualizeDeformedMesh(obj, uGlobal, scale, titleStr)
            crd = obj.globalMesh.coord;
            con = obj.globalMesh.connec;
            uResh = reshape(uGlobal, obj.globalMesh.ndim, size(crd, 1))';
            defCrd = crd + scale * uResh;
            dispMag = sqrt(uResh(:, 1).^2 + uResh(:, 2).^2);
            
            figure('Name', titleStr, 'Color', 'w'); hold on; axis equal;
            patch('Faces', con, 'Vertices', crd, 'FaceColor', 'none', 'EdgeColor', [0.8 0.8 0.8], 'LineStyle', '--');
            patch('Faces', con, 'Vertices', defCrd, 'FaceVertexCData', dispMag, 'FaceColor', 'interp', 'EdgeColor', '#333333');
            colormap(jet); colorbar; title(titleStr);
        end
    end
end