classdef tempFetiDPElasticityNofMat < handle
    
    properties (Access = private)
        localStiffness
        localForces
        localMeshes
        meshCoords
        dofsPerNode
        nodeTol
        numSubdomains
        boundaryConditions
        
        primalDofsGlobal
        dualDofsGlobal
        remDofsGlobal
        
        primalDofsLocal
        dualDofsLocal
        remDofsLocal
        
        primalIdxLocal
        dualIdxLocal
        dualSignsLocal
        
        krrFactors
        tdrLocal
        sppActive
        brActive
        kiiFactors
    end
    
    methods (Access = public)
        
        % -----------------------------------------------------------------
        % 1. CONSTRUCTOR
        % -----------------------------------------------------------------
        function obj = tempFetiDPElasticityNofMat(globalMesh, subMeshes, stiffness, forces, tol, dofsNode, boundaryConditions, localToGlobalMaps)
            obj.meshCoords         = globalMesh.coord;
            obj.dofsPerNode        = dofsNode;
            obj.localMeshes        = subMeshes;
            obj.localStiffness     = stiffness;
            obj.localForces        = forces;
            obj.nodeTol            = tol;
            obj.numSubdomains      = length(stiffness);
            obj.boundaryConditions = boundaryConditions;
            obj.primalIdxLocal     = cell(obj.numSubdomains, 1);
            
            obj.extractFetiDofs(localToGlobalMaps);
        end
        
        % -----------------------------------------------------------------
        % 2. MAIN PROBLEM ASSEMBLY
        % -----------------------------------------------------------------
        function [fMatOperator, dBar] = assembleProblem(obj)
            allPrimals = unique(vertcat(obj.primalDofsGlobal{:}));
            allDuals   = setdiff(unique(vertcat(obj.dualDofsGlobal{:})), allPrimals);
            
            numPrimals = length(allPrimals);
            numDuals   = length(allDuals);
            
            dBar        = zeros(numDuals, 1);
            sppBase     = sparse(numPrimals, numPrimals);
            brKrrInvKrp = sparse(numDuals, numPrimals);
            rhsPrimal   = zeros(numPrimals, 1);
            
            visitedDuals = zeros(max(allDuals), 1);
            
            obj.krrFactors = cell(obj.numSubdomains, 1);
            obj.kiiFactors = cell(obj.numSubdomains, 1);
            obj.tdrLocal   = cell(obj.numSubdomains, 1);
            
            for subId = 1:obj.numSubdomains
                [Krr, Krp, Kpr, Kpp, fR, fP, Tdr] = obj.splitLocalMatrices(subId);
                
                % Scalability core: Precompute LU/Cholesky factorizations
                decKrr = decomposition(Krr);
                obj.krrFactors{subId} = decKrr;
                obj.tdrLocal{subId}   = Tdr;
                
                rLoc = obj.remDofsLocal{subId};
                nRem = length(rLoc);
                obj.kiiFactors{subId} = decomposition(Krr(1:nRem, 1:nRem));
                
                Urd = decKrr \ Tdr';
                Urp = decKrr \ Krp;
                
                pGlobal = obj.primalDofsGlobal{subId};
                [~, pRows] = ismember(pGlobal, allPrimals);
                obj.primalIdxLocal{subId} = pRows;
                
                dGlobal       = obj.dualDofsGlobal{subId};
                [~, dualRows] = ismember(dGlobal, allDuals);
                
                isFirst               = ~visitedDuals(dGlobal);
                dualSigns             = ones(length(dGlobal), 1);
                dualSigns(~isFirst)   = -1;
                visitedDuals(dGlobal) = true;
                
                obj.dualIdxLocal{subId}   = dualRows;
                obj.dualSignsLocal{subId} = dualSigns;
                
                sppLoc = Kpp - Kpr * Urp;
                sppBase(pRows, pRows) = sppBase(pRows, pRows) + sppLoc;
                
                if ~isempty(dGlobal)
                    dBar(dualRows)               = dBar(dualRows) + dualSigns .* (Urd' * fR);
                    brKrrInvKrp(dualRows, pRows) = brKrrInvKrp(dualRows, pRows) + dualSigns .* (Tdr * Urp);
                end
                
                rhsPrimal(pRows) = rhsPrimal(pRows) + (fP - Urp' * fR);
            end
            
            activeDofs    = obj.getActivePrimalDofs(allPrimals);
            obj.sppActive = sppBase(activeDofs, activeDofs);
            obj.brActive  = brKrrInvKrp(:, activeDofs);
            rhsActive     = rhsPrimal(activeDofs);
            
            dBar = dBar - obj.brActive * (obj.sppActive \ rhsActive);
            
            % Matrix-free operator return
            fMatOperator = @(lambda) obj.applyGlobalF(lambda);
        end
        
        function fOut = applyGlobalF(obj, lambda)
            numDuals = length(lambda);
            fOut     = zeros(numDuals, 1);
            
            for subId = 1:obj.numSubdomains
                dualRows  = obj.dualIdxLocal{subId};
                dualSigns = obj.dualSignsLocal{subId};
                
                if isempty(dualRows)
                    continue;
                end
                
                lambdaLoc = dualSigns .* lambda(dualRows);
                Tdr       = obj.tdrLocal{subId};
                decKrr    = obj.krrFactors{subId};
                
                termLoc = Tdr * (decKrr \ (Tdr' * lambdaLoc));
                fOut(dualRows) = fOut(dualRows) + dualSigns .* termLoc;
            end
            
            primalCorrection = obj.brActive * (obj.sppActive \ (obj.brActive' * lambda));
            fOut = fOut + primalCorrection;
        end
        
        % -----------------------------------------------------------------
        % 3. SOLVER TOOLS (PRECONDITIONER)
        % -----------------------------------------------------------------
        function M = buildPrecondMatrix(obj)
            allPrimals = unique(vertcat(obj.primalDofsGlobal{:}));
            allDuals   = setdiff(unique(vertcat(obj.dualDofsGlobal{:})), allPrimals);
            numDuals   = length(allDuals);
            
            multiplicity = zeros(numDuals, 1);
            for subId = 1:obj.numSubdomains
                dualRows = obj.dualIdxLocal{subId};
                multiplicity(dualRows) = multiplicity(dualRows) + 1;
            end
            
            M = zeros(numDuals);
            for subId = 1:obj.numSubdomains
                dualRows  = obj.dualIdxLocal{subId};
                dualSigns = obj.dualSignsLocal{subId};
                if isempty(dualRows)
                    continue;
                end
                
                Sdd = obj.computeLocalSchur(subId);
                w   = 1 ./ multiplicity(dualRows);
                
                wSigned = diag(w .* dualSigns);
                M(dualRows, dualRows) = M(dualRows, dualRows) + wSigned * Sdd * wSigned;
            end
        end
        
        function z = applyDirichletPrecond(obj, r)
            allPrimals = unique(vertcat(obj.primalDofsGlobal{:}));
            allDuals   = setdiff(unique(vertcat(obj.dualDofsGlobal{:})), allPrimals);
            numDuals   = length(allDuals);
            
            multiplicity = zeros(numDuals, 1);
            for subId = 1:obj.numSubdomains
                dualRows = obj.dualIdxLocal{subId};
                multiplicity(dualRows) = multiplicity(dualRows) + 1;
            end
            
            z = zeros(numDuals, 1);
            for subId = 1:obj.numSubdomains
                dualRows  = obj.dualIdxLocal{subId};
                dualSigns = obj.dualSignsLocal{subId};
                if isempty(dualRows)
                    continue;
                end
                
                w    = 1 ./ multiplicity(dualRows);
                rLoc = dualSigns .* r(dualRows);
                
                sddLoc = obj.applyLocalSchurVector(subId, w .* rLoc);
                z(dualRows) = z(dualRows) + dualSigns .* (w .* sddLoc);
            end
        end
        
        % -----------------------------------------------------------------
        % 4. POST-PROCESSING
        % -----------------------------------------------------------------
        function uFull = reconstructGlobalSolution(obj, lambdaSol, numNodes)
            allPrimals = unique(vertcat(obj.primalDofsGlobal{:}));
            numP       = length(allPrimals);
            sppBase    = sparse(numP, numP);
            rhsPrimal  = zeros(numP, 1);
            
            urpCell      = cell(obj.numSubdomains, 1);
            urdCell      = cell(obj.numSubdomains, 1);
            krrInvFrCell = cell(obj.numSubdomains, 1);
            
            for subId = 1:obj.numSubdomains
                [~, Krp, Kpr, Kpp, fR, fP, Tdr] = obj.splitLocalMatrices(subId);
                pRows = obj.getPrimalRows(subId, allPrimals);
                
                decKrr = obj.krrFactors{subId};
                
                Urp      = decKrr \ Krp;
                Urd      = decKrr \ Tdr';
                KrrInvFr = decKrr \ fR;
                
                urpCell{subId}      = Urp;
                urdCell{subId}      = Urd;
                krrInvFrCell{subId} = KrrInvFr;
                
                sppLoc = Kpp - Kpr * Urp;
                sppBase(pRows, pRows) = sppBase(pRows, pRows) + sppLoc;
                
                dualRows  = obj.dualIdxLocal{subId};
                dualSigns = obj.dualSignsLocal{subId};
                lambdaLoc = dualSigns .* lambdaSol(dualRows);
                term      = Tdr' * lambdaLoc - fR;
                rhsPrimal(pRows) = rhsPrimal(pRows) + (fP + Urp' * term);
            end
            
            activeDofs = obj.getActivePrimalDofs(allPrimals);
            sppActive  = sppBase(activeDofs, activeDofs);
            rhsActive  = rhsPrimal(activeDofs);
            
            uP = zeros(numP, 1);
            uP(activeDofs) = sppActive \ rhsActive;
            
            uFull = zeros(numNodes * obj.dofsPerNode, 1);
            
            for subId = 1:obj.numSubdomains
                pGlobal = obj.primalDofsGlobal{subId};
                dGlobal = obj.dualDofsGlobal{subId};
                rGlobal = obj.remDofsGlobal{subId};
                
                pLocal = obj.primalDofsLocal{subId};
                dLocal = obj.dualDofsLocal{subId};
                rLocal = obj.remDofsLocal{subId};
                
                pRows = obj.getPrimalRows(subId, allPrimals);
                
                dualRows  = obj.dualIdxLocal{subId};
                dualSigns = obj.dualSignsLocal{subId};
                lambdaLoc = dualSigns .* lambdaSol(dualRows);
                
                uPLoc   = uP(pRows);
                uRemLoc = krrInvFrCell{subId} - urpCell{subId} * uPLoc - urdCell{subId} * lambdaLoc;
                
                numLocalDofs = length(pLocal) + length(dLocal) + length(rLocal);
                uLocalTilde = zeros(numLocalDofs, 1);
                uLocalTilde(pLocal) = uPLoc;
                uLocalTilde([rLocal; dLocal]) = uRemLoc;
                
                uFull(pGlobal)            = uLocalTilde(pLocal);
                uFull([rGlobal; dGlobal]) = uLocalTilde([rLocal; dLocal]);
            end
        end
        
        % -----------------------------------------------------------------
        % 5. VISUALIZATION
        % -----------------------------------------------------------------
        function visualizeFetiNodes(obj)
            allPrimalDofs = unique(vertcat(obj.primalDofsGlobal{:}));
            allDualDofs   = unique(vertcat(obj.dualDofsGlobal{:}));
            allRemDofs    = unique(vertcat(obj.remDofsGlobal{:}));
            
            primalNodes = unique(ceil(allPrimalDofs / obj.dofsPerNode));
            dualNodes   = unique(ceil(allDualDofs   / obj.dofsPerNode));
            remNodes    = unique(ceil(allRemDofs    / obj.dofsPerNode));
            
            dualNodes = setdiff(dualNodes, primalNodes);
            remNodes  = setdiff(remNodes, union(primalNodes, dualNodes));
            
            pCoords = obj.meshCoords(primalNodes, :);
            dCoords = obj.meshCoords(dualNodes, :);
            rCoords = obj.meshCoords(remNodes, :);
            
            figure('Name', 'FETI-DP Nodes (Elasticity)', 'Color', 'w');
            hold on; axis equal;
            
            for i = 1:obj.numSubdomains
                patch('Faces', obj.localMeshes{i}.connec, 'Vertices', obj.localMeshes{i}.coord, ...
                    'FaceColor', 'none', 'EdgeColor', [0.8 0.8 0.8], 'LineWidth', 1.5, 'HandleVisibility', 'off');
            end
            
            scatter(rCoords(:, 1), rCoords(:, 2), 20, [0.5 0.5 0.5], 'filled', 'DisplayName', 'Interior');
            scatter(dCoords(:, 1), dCoords(:, 2), 40, 'b', 'filled', 'DisplayName', 'Interface (Dual)');
            scatter(pCoords(:, 1), pCoords(:, 2), 80, 'r', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1, 'DisplayName', 'Corners/Edges (Primal)');
            
            legend('Location', 'bestoutside');
            title('FETI-DP Node Distribution (2D Elasticity)');
            xlabel('X'); ylabel('Y');
            grid on; hold off;
        end
    end
    
    methods (Access = private)
        
        % -----------------------------------------------------------------
        % 6. TOPOLOGICAL SETUP
        % -----------------------------------------------------------------
        function extractFetiDofs(obj, localToGlobalMaps)
            nSub = obj.numSubdomains;

            pGlobal = cell(nSub, 1); 
            dGlobal = cell(nSub, 1); 
            rGlobal = cell(nSub, 1);
            pLocal  = cell(nSub, 1); 
            dLocal  = cell(nSub, 1); 
            rLocal  = cell(nSub, 1);
            
            numGlobalNodes   = size(obj.meshCoords, 1);

            % nodeSubMask      = sparse(numGlobalNodes, nSub);
            % globalNodesCache = cell(nSub, 1);
            % 
            % for i = 1:nSub
            %     gN = localToGlobalMaps{i};
            %     nodeSubMask(gN, i) = 1;
            % end
            nodeMultiplicity = zeros(numGlobalNodes, 1);
            for i = 1:nSub
                gN = localToGlobalMaps{i};
                nodeMultiplicity(gN) = nodeMultiplicity(gN) + 1;
            end

            for i = 1:nSub
                gNodes = localToGlobalMaps{i};
                nNodes = length(gNodes);
                
                % Primal Nodes
                boundaryMeshes = obj.localMeshes{i}.createBoundaryMesh();
                numBoundaryMeshes = length(boundaryMeshes);
                primalNodesPerMesh = cell(numBoundaryMeshes, 1);
                
                for b = 1:numBoundaryMeshes
                    boundaryConnectivity = boundaryMeshes{b}.globalConnec(:);                    
                    nodeOccurrences = accumarray(boundaryConnectivity, 1);
                    primalNodesPerMesh{b} = find(nodeOccurrences == 1);
                end                
                idxPrimal = unique(vertcat(primalNodesPerMesh{:}));
                
                % Dual
                isShared = nodeMultiplicity(gNodes) > 1;
                idxDualFinal = setdiff(find(isShared), idxPrimal);
                
                idxRem = setdiff((1:nNodes)', [idxPrimal; idxDualFinal]);
                
                pGlobal{i} = obj.nodesToDofs(gNodes(idxPrimal));
                dGlobal{i} = obj.nodesToDofs(gNodes(idxDualFinal));
                rGlobal{i} = obj.nodesToDofs(gNodes(idxRem));
                
                pLocal{i} = obj.nodesToDofs(idxPrimal);
                dLocal{i} = obj.nodesToDofs(idxDualFinal);
                rLocal{i} = obj.nodesToDofs(idxRem);
            end
            
            obj.primalDofsGlobal = pGlobal; 
            obj.primalDofsLocal  = pLocal;
            obj.dualDofsGlobal   = dGlobal; 
            obj.dualDofsLocal    = dLocal;
            obj.remDofsGlobal    = rGlobal; 
            obj.remDofsLocal     = rLocal;
        end
        
        function dofs = nodesToDofs(obj, nodes)
            nodes = nodes(:);
            dofs  = zeros(length(nodes) * obj.dofsPerNode, 1);
            for d = 1:obj.dofsPerNode
                dofs(d:obj.dofsPerNode:end) = (nodes - 1) * obj.dofsPerNode + d;
            end
        end
        
        % -----------------------------------------------------------------
        % 7. ALGEBRAIC TRANSFORMATIONS
        % -----------------------------------------------------------------
        function [Krr, Krp, Kpr, Kpp, fR, fP, Tdr] = splitLocalMatrices(obj, subId)
            kMat = obj.localStiffness{subId};
            fVec = obj.localForces{subId};
            
            pLoc = obj.primalDofsLocal{subId};
            dLoc = obj.dualDofsLocal{subId};
            rLoc = obj.remDofsLocal{subId};
            
            remLoc = [rLoc; dLoc];
            
            Krr = kMat(remLoc, remLoc);
            Krp = kMat(remLoc, pLoc);
            Kpr = kMat(pLoc, remLoc);
            Kpp = kMat(pLoc, pLoc);
            
            fR  = fVec(remLoc);
            fP  = fVec(pLoc);
            
            numD = length(dLoc);
            numR = length(remLoc);
            
            rows = (1:numD)';
            cols = (length(rLoc) + 1 : numR)';
            vals = ones(numD, 1);
            Tdr  = sparse(rows, cols, vals, numD, numR);
        end
        
        % -----------------------------------------------------------------
        % 8. MATHEMATICAL OPERATIONS
        % -----------------------------------------------------------------
        function sddX = applyLocalSchurVector(obj, subId, x)
            kMat = obj.localStiffness{subId};
            rLoc = obj.remDofsLocal{subId};
            dLoc = obj.dualDofsLocal{subId};
            
            Kid = kMat(rLoc, dLoc);
            Kdi = kMat(dLoc, rLoc);
            Kdd = kMat(dLoc, dLoc);
            
            decKii = obj.kiiFactors{subId};
            sddX   = Kdd * x - Kdi * (decKii \ (Kid * x));
        end
        
        function Sdd = computeLocalSchur(obj, subId)
            kMat = obj.localStiffness{subId};
            rLoc = obj.remDofsLocal{subId};
            dLoc = obj.dualDofsLocal{subId};
            
            Kii = kMat(rLoc, rLoc);
            Kid = kMat(rLoc, dLoc);
            Kdi = kMat(dLoc, rLoc);
            Kdd = kMat(dLoc, dLoc);
            
            Sdd = Kdd - Kdi * (Kii \ Kid);
        end
        
        % -----------------------------------------------------------------
        % 9. GETTERS & INDEXING
        % -----------------------------------------------------------------
        function pRows = getPrimalRows(obj, subId, allPrimals)
            pRows = obj.primalIdxLocal{subId};
            
            if isempty(pRows)
                [~, pRows] = ismember(obj.primalDofsGlobal{subId}, allPrimals);
            end
        end
        
        function activeIdx = getActivePrimalDofs(obj, allPrimals)
            activeIdx = (1:length(allPrimals));
        end
    end
end