classdef tempFetiDPElasticity3DEdgeAverage < handle
    
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
        edgeDofsGrouped
        
        B_matrices
        W_cell
        numMultipliers
    end
    
    methods (Access = public)
        
        % -----------------------------------------------------------------
        % 1. CONSTRUCTOR
        % -----------------------------------------------------------------
        function obj = tempFetiDPElasticity3DEdgeAverage(globalMesh, subMeshes, stiffness, forces, tol, dofsNode, boundaryConditions)
            obj.meshCoords         = globalMesh.coord;
            obj.dofsPerNode        = dofsNode;
            obj.localMeshes        = subMeshes;
            obj.localStiffness     = stiffness;
            obj.localForces        = forces;
            obj.nodeTol            = tol;
            obj.numSubdomains      = length(stiffness);
            obj.boundaryConditions = boundaryConditions;
            obj.primalIdxLocal     = cell(obj.numSubdomains, 1);
            obj.edgeDofsGrouped    = cell(obj.numSubdomains, 1);
            
            obj.extractFetiDofs();
            obj.buildBooleanMatrices();
        end
        
        % -----------------------------------------------------------------
        % 2. MAIN PROBLEM ASSEMBLY
        % -----------------------------------------------------------------
        function [fMat, dBar] = assembleProblem(obj)
            allPrimals = unique(vertcat(obj.primalDofsGlobal{:}));
            numPrimals = length(allPrimals);
            numDuals   = obj.numMultipliers;
            
            fDual       = sparse(numDuals, numDuals);
            dBar        = zeros(numDuals, 1);
            SPP         = sparse(numPrimals, numPrimals);
            BrKrrInvKrp = sparse(numDuals, numPrimals);
            rhsPrimal   = zeros(numPrimals, 1);
            
            for subId = 1:obj.numSubdomains
                [Krr, Krp, Kpr, Kpp, fR, fP, Tdr] = obj.splitLocalMatrices(subId);
                
                Urd = Krr \ Tdr';
                Urp = Krr \ Krp;
                
                pGlobal = obj.primalDofsGlobal{subId};
                [~, pRows] = ismember(pGlobal, allPrimals);
                obj.primalIdxLocal{subId} = pRows;
                
                sppLoc = Kpp - Kpr * Urp;
                SPP(pRows, pRows) = SPP(pRows, pRows) + sppLoc;
                
                B_sub = obj.B_matrices{subId};
                if nnz(B_sub) > 0
                    F_local = Tdr * Urd;
                    fDual = fDual + B_sub * F_local * B_sub';
                    dBar = dBar + B_sub * (Tdr * (Krr \ fR));
                    BrKrrInvKrp(:, pRows) = BrKrrInvKrp(:, pRows) + B_sub * (Tdr * Urp);
                end
                
                rhsPrimal(pRows) = rhsPrimal(pRows) + (fP - Urp' * fR);
            end
            
            activeDofs = obj.getActivePrimalDofs(allPrimals);
            SPPActive  = SPP(activeDofs, activeDofs);
            BrActive   = BrKrrInvKrp(:, activeDofs);
            rhsActive  = rhsPrimal(activeDofs);
            
            fMat = fDual + BrActive * (SPPActive \ BrActive');
            dBar = dBar  - BrActive * (SPPActive \ rhsActive);
        end
        
        % -----------------------------------------------------------------
        % 3. SOLVER TOOLS (PRECONDITIONER)
        % -----------------------------------------------------------------
        function z = applyDirichletPrecond(obj, r)
            z = zeros(obj.numMultipliers, 1);
            for subId = 1:obj.numSubdomains
                B_sub = obj.B_matrices{subId};
                if nnz(B_sub) == 0, continue; end
                
                Sdd = obj.computeLocalSchur(subId);
                W   = diag(obj.W_cell{subId});
                
                rLoc = B_sub' * r;
                zLoc = W * Sdd * W * rLoc;
                z = z + B_sub * zLoc;
            end
        end
        
        function M = buildPrecondMatrix(obj)
            M = zeros(obj.numMultipliers, obj.numMultipliers);
            for subId = 1:obj.numSubdomains
                B_sub = obj.B_matrices{subId};
                if nnz(B_sub) == 0, continue; end
                
                Sdd = obj.computeLocalSchur(subId);
                W   = diag(obj.W_cell{subId});
                
                M = M + B_sub * W * Sdd * W * B_sub';
            end
        end
        
        % -----------------------------------------------------------------
        % 4. POST-PROCESSING
        % -----------------------------------------------------------------
        function uFull = reconstructGlobalSolution(obj, lambdaSol, numNodes)
            allPrimals = unique(vertcat(obj.primalDofsGlobal{:}));
            numP       = length(allPrimals);
            SPP        = sparse(numP, numP);
            rhsPrimal  = zeros(numP, 1);
            
            UrpCell      = cell(obj.numSubdomains, 1);
            UrdCell      = cell(obj.numSubdomains, 1);
            KrrInvFrCell = cell(obj.numSubdomains, 1);
            
            for subId = 1:obj.numSubdomains
                [Krr, Krp, Kpr, Kpp, fR, fP, Tdr] = obj.splitLocalMatrices(subId);
                pRows = obj.getPrimalRows(subId, allPrimals);
                
                Urp      = Krr \ Krp;
                Urd      = Krr \ Tdr';
                KrrInvFr = Krr \ fR;
                
                UrpCell{subId}      = Urp;
                UrdCell{subId}      = Urd;
                KrrInvFrCell{subId} = KrrInvFr;
                
                sppLoc = Kpp - Kpr * Urp;
                SPP(pRows, pRows) = SPP(pRows, pRows) + sppLoc;
                
                B_sub = obj.B_matrices{subId};
                if nnz(B_sub) > 0
                    lambdaLoc = B_sub' * lambdaSol;
                else
                    lambdaLoc = zeros(length(obj.dualDofsGlobal{subId}), 1);
                end
                
                term = Tdr' * lambdaLoc - fR;
                rhsPrimal(pRows) = rhsPrimal(pRows) + (fP + Urp' * term);
            end
            
            activeDofs = obj.getActivePrimalDofs(allPrimals);
            uP = zeros(numP, 1);
            uP(activeDofs) = SPP(activeDofs, activeDofs) \ rhsPrimal(activeDofs);
            
            uFull = zeros(numNodes * obj.dofsPerNode, 1);
            
            for subId = 1:obj.numSubdomains
                pGlobal = obj.primalDofsGlobal{subId};
                dGlobal = obj.dualDofsGlobal{subId};
                rGlobal = obj.remDofsGlobal{subId};
                
                pLocal = obj.primalDofsLocal{subId};
                dLocal = obj.dualDofsLocal{subId};
                rLocal = obj.remDofsLocal{subId};
                
                pRows = obj.getPrimalRows(subId, allPrimals);
                
                B_sub = obj.B_matrices{subId};
                if nnz(B_sub) > 0
                    lambdaLoc = B_sub' * lambdaSol;
                else
                    lambdaLoc = zeros(length(dLocal), 1);
                end
                
                uPLoc   = uP(pRows);
                uRemLoc = KrrInvFrCell{subId} - UrpCell{subId} * uPLoc - UrdCell{subId} * lambdaLoc;
                
                numLocalDofs = length(pLocal) + length(dLocal) + length(rLocal);
                uLocalTilde = zeros(numLocalDofs, 1);
                uLocalTilde(pLocal) = uPLoc;
                uLocalTilde([rLocal; dLocal]) = uRemLoc;
                
                edgeGroups = obj.edgeDofsGrouped{subId};
                uLocalPhys = obj.applyEdgeAverageForward(uLocalTilde, edgeGroups);
                
                uFull(pGlobal)            = uLocalPhys(pLocal);
                uFull([rGlobal; dGlobal]) = uLocalPhys([rLocal; dLocal]);
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
            
            figure('Name', 'FETI-DP Nodes (3D Elasticity)', 'Color', 'w');
            hold on; axis equal;
            
            scatter3(obj.meshCoords(remNodes, 1), obj.meshCoords(remNodes, 2), obj.meshCoords(remNodes, 3), ...
                16, [0.5 0.5 0.5], 'filled', 'DisplayName', 'Interior');
            scatter3(obj.meshCoords(dualNodes, 1), obj.meshCoords(dualNodes, 2), obj.meshCoords(dualNodes, 3), ...
                32, 'b', 'filled', 'DisplayName', 'Interface (Dual)');
            scatter3(obj.meshCoords(primalNodes, 1), obj.meshCoords(primalNodes, 2), obj.meshCoords(primalNodes, 3), ...
                70, 'r', 'filled', 'MarkerEdgeColor', 'k', 'DisplayName', 'Corners/Edges (Primal)');
            
            legend('Location', 'bestoutside');
            title('FETI-DP Node Distribution (3D)');
            xlabel('X'); ylabel('Y'); zlabel('Z');
            grid on; view(3); hold off;
        end
    end
    
    methods (Access = private)
        
        % -----------------------------------------------------------------
        % 6. TOPOLOGICAL SETUP
        % -----------------------------------------------------------------
        function extractFetiDofs(obj)
            nSub = obj.numSubdomains;
            tol  = obj.nodeTol;
            
            pGlobal = cell(nSub, 1); 
            dGlobal = cell(nSub, 1); 
            rGlobal = cell(nSub, 1);
            pLocal  = cell(nSub, 1); 
            dLocal  = cell(nSub, 1); 
            rLocal  = cell(nSub, 1);
            
            numGlobalNodes   = size(obj.meshCoords, 1);
            nodeSubMask      = sparse(numGlobalNodes, nSub);
            globalNodesCache = cell(nSub, 1);
            
            for i = 1:nSub
                [~, gN] = ismembertol(obj.localMeshes{i}.coord, obj.meshCoords, tol, 'ByRows', true);
                globalNodesCache{i} = gN;
                nodeSubMask(gN, i) = 1;
            end
            
            for i = 1:nSub
                coords = obj.localMeshes{i}.coord;
                gNodes = globalNodesCache{i};
                nNodes = size(coords, 1);
                
                idxCorner    = obj.computeLocalCornerNodes(coords);
                idxPrimal    = idxCorner;
                idxDualFinal = [];
                
                edgeNodesByLine  = obj.computeLocalEdgeNodeSets(coords, idxCorner);
                numEdges         = length(edgeNodesByLine);
                edgeDofsLocalSub = cell(numEdges, 1);
                edgeNodeSet      = [];
                
                for e = 1:numEdges
                    localNodesInEdge = edgeNodesByLine{e};
                    if isempty(localNodesInEdge)
                        continue;
                    end
                    
                    % Retain only edge nodes that are genuinely shared
                    edgeGlobalNodes  = gNodes(localNodesInEdge);
                    isShared         = sum(nodeSubMask(edgeGlobalNodes, :), 2) > 1;
                    localNodesInEdge = localNodesInEdge(isShared);
                    
                    if isempty(localNodesInEdge)
                        continue;
                    end
                    
                    edgeNodeSet = [edgeNodeSet; localNodesInEdge(:)];
                    
                    [~, order]       = sort(gNodes(localNodesInEdge));
                    localNodesInEdge = localNodesInEdge(order);
                    
                    if length(localNodesInEdge) > 1
                        eDofs = obj.nodesToDofs(localNodesInEdge);
                        edgeDofsLocalSub{e} = eDofs;
                        
                        idxPrimal    = [idxPrimal; localNodesInEdge(1)];
                        idxDualFinal = [idxDualFinal; localNodesInEdge(2:end)];
                    else
                        idxDualFinal = [idxDualFinal; localNodesInEdge(1)];
                    end
                end
                
                sharedNodes   = find(sum(nodeSubMask(gNodes, :), 2) > 1);
                faceDualNodes = setdiff(sharedNodes, unique([idxCorner; edgeNodeSet]));
                idxDualFinal  = [idxDualFinal; faceDualNodes];
                
                idxRem = setdiff((1:nNodes)', [idxPrimal; idxDualFinal]);
                
                pGlobal{i} = obj.nodesToDofs(gNodes(idxPrimal));
                dGlobal{i} = obj.nodesToDofs(gNodes(idxDualFinal));
                rGlobal{i} = obj.nodesToDofs(gNodes(idxRem));
                
                pLocal{i} = obj.nodesToDofs(idxPrimal);
                dLocal{i} = obj.nodesToDofs(idxDualFinal);
                rLocal{i} = obj.nodesToDofs(idxRem);
                
                obj.edgeDofsGrouped{i} = edgeDofsLocalSub(~cellfun('isempty', edgeDofsLocalSub));
            end
            
            obj.primalDofsGlobal = pGlobal; obj.primalDofsLocal = pLocal;
            obj.dualDofsGlobal   = dGlobal; obj.dualDofsLocal   = dLocal;
            obj.remDofsGlobal    = rGlobal; obj.remDofsLocal    = rLocal;
        end
        
        function buildBooleanMatrices(obj)
            allDuals = [];
            subList = [];
            localDualIdxList = [];
            
            for subId = 1:obj.numSubdomains
                dG = obj.dualDofsGlobal{subId};
                if isempty(dG)
                    continue; 
                end
                allDuals = [allDuals; dG];
                subList = [subList; subId * ones(length(dG), 1)];
                localDualIdxList = [localDualIdxList; (1:length(dG))'];
            end
            
            if isempty(allDuals)
                obj.numMultipliers = 0;
                obj.B_matrices = cell(obj.numSubdomains, 1);
                obj.W_cell = cell(obj.numSubdomains, 1);
                for s = 1:obj.numSubdomains
                    obj.B_matrices{s} = sparse(0, length(obj.dualDofsGlobal{s}));
                    obj.W_cell{s} = [];
                end
                return;
            end
            
            [sortedDuals, sortIdx] = sort(allDuals);
            sortedSubs = subList(sortIdx);
            sortedLocalIdx = localDualIdxList(sortIdx);
            
            [uniqueDuals, firstIdx, ~] = unique(sortedDuals, 'first');
            [~, lastIdx, ~] = unique(sortedDuals, 'last');
            
            obj.numMultipliers = 0;
            B_rows = cell(obj.numSubdomains, 1);
            B_cols = cell(obj.numSubdomains, 1);
            B_vals = cell(obj.numSubdomains, 1);
            obj.W_cell = cell(obj.numSubdomains, 1);
            
            for s = 1:obj.numSubdomains
                nLoc = length(obj.dualDofsGlobal{s});
                obj.W_cell{s} = zeros(nLoc, 1);
            end
            
            for i = 1:length(uniqueDuals)
                idxRange = firstIdx(i):lastIdx(i);
                subs = sortedSubs(idxRange);
                locs = sortedLocalIdx(idxRange);
                m = length(subs);
                
                % Multiplicidad para el precondicionador
                for j = 1:m
                    obj.W_cell{subs(j)}(locs(j)) = 1 / m;
                end
                
                % Multiplicadores: Unimos el primero de la lista con todos los demas
                s1 = subs(1);
                loc1 = locs(1);
                for j = 2:m
                    sj = subs(j);
                    locj = locs(j);
                    obj.numMultipliers = obj.numMultipliers + 1;
                    
                    B_rows{s1} = [B_rows{s1}; obj.numMultipliers];
                    B_cols{s1} = [B_cols{s1}; loc1];
                    B_vals{s1} = [B_vals{s1}; 1];
                    
                    B_rows{sj} = [B_rows{sj}; obj.numMultipliers];
                    B_cols{sj} = [B_cols{sj}; locj];
                    B_vals{sj} = [B_vals{sj}; -1];
                end
            end
            
            obj.B_matrices = cell(obj.numSubdomains, 1);
            for s = 1:obj.numSubdomains
                nLoc = length(obj.dualDofsGlobal{s});
                if isempty(B_rows{s})
                    obj.B_matrices{s} = sparse(obj.numMultipliers, nLoc);
                else
                    obj.B_matrices{s} = sparse(B_rows{s}, B_cols{s}, B_vals{s}, obj.numMultipliers, nLoc);
                end
            end
        end
        
        function idxCorner = computeLocalCornerNodes(obj, coords)
            mins = min(coords, [], 1);
            maxs = max(coords, [], 1);
            isExtreme = abs(coords - mins) < obj.nodeTol | abs(coords - maxs) < obj.nodeTol;
            idxCorner = find(all(isExtreme, 2));
        end

        function edgeSets = computeLocalEdgeNodeSets(obj, coords, idxCorner)
            mins = min(coords, [], 1);
            maxs = max(coords, [], 1);
            edgeSets = {};
            iEdge = 0;
            for freeDim = 1:3
                fixedDims = setdiff(1:3, freeDim);
                for s1 = 1:2
                    for s2 = 1:2
                        vals = [mins(fixedDims(1)), mins(fixedDims(2))];
                        if s1 == 2
                            vals(1) = maxs(fixedDims(1));
                        end
                        if s2 == 2
                            vals(2) = maxs(fixedDims(2));
                        end
                        isOnEdge = abs(coords(:, fixedDims(1)) - vals(1)) < obj.nodeTol & ...
                                   abs(coords(:, fixedDims(2)) - vals(2)) < obj.nodeTol;
                        edgeNodes = setdiff(find(isOnEdge), idxCorner);
                        [~, order] = sort(coords(edgeNodes, freeDim));
                        
                        iEdge = iEdge + 1;
                        edgeSets{iEdge} = edgeNodes(order);
                    end
                end
            end
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
            
            edgeGroups   = obj.edgeDofsGrouped{subId};
            [kMat, fVec] = obj.applyEdgeAverageTransformation(kMat, fVec, edgeGroups);
            
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
            if numD == 0
                rows = []; cols = []; vals = [];
            end
            Tdr  = sparse(rows, cols, vals, numD, numR);
        end
        
        function [KMod, fMod] = applyEdgeAverageTransformation(obj, K, f, edgeGroups)
            numDofs = size(K, 1);
            T = speye(numDofs);
            
            for e = 1:length(edgeGroups)
                eDofs = edgeGroups{e};
                for d = 1:obj.dofsPerNode
                    dimDofs = eDofs(d:obj.dofsPerNode:end);
                    pAvg    = dimDofs(1);
                    dDeltas = dimDofs(2:end);
                    
                    T(dimDofs, pAvg) = 1;
                    T(pAvg, dDeltas) = -1;
                end
            end
            
            KMod = T' * K * T;
            fMod = T' * f;
        end
        
        function uPhys = applyEdgeAverageForward(obj, uTilde, edgeGroups)
            numDofs = length(uTilde);
            T = speye(numDofs);
            
            for e = 1:length(edgeGroups)
                eDofs = edgeGroups{e};
                for d = 1:obj.dofsPerNode
                    dimDofs = eDofs(d:obj.dofsPerNode:end);
                    pAvg    = dimDofs(1);
                    dDeltas = dimDofs(2:end);
                    
                    T(dimDofs, pAvg) = 1;
                    T(pAvg, dDeltas) = -1;
                end
            end
            
            uPhys = T * uTilde;
        end
        
        % -----------------------------------------------------------------
        % 8. MATHEMATICAL OPERATIONS
        % -----------------------------------------------------------------
        function Sdd = computeLocalSchur(obj, subId)
            kMat = obj.localStiffness{subId};
            rLoc = obj.remDofsLocal{subId};
            dLoc = obj.dualDofsLocal{subId};
            
            edgeGroups = obj.edgeDofsGrouped{subId};
            fVecDummy  = zeros(size(kMat, 1), 1);
            [kMat, ~]  = obj.applyEdgeAverageTransformation(kMat, fVecDummy, edgeGroups);
            
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
        
        function activeIdx = getActivePrimalDofs(~, allPrimals)
            activeIdx = (1:length(allPrimals));
        end
    end
end