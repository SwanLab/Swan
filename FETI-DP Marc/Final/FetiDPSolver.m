classdef FetiDPSolver < handle
    % FETI-DP algebraic solver for 2D general physics on a domain decomposition.    

    properties (Access = private)
        localStiffness
        localForces
        localMeshes
        meshCoords
        dofsPerNode
        nodeTol
        numSubdomains
        boundaryConditions
        
        useMatrixFree
        useEdgeAverage
        
        primalDofsGlobal
        dualDofsGlobal
        remDofsGlobal
        primalDofsLocal
        dualDofsLocal
        remDofsLocal
        
        primalIdxLocal
        dualIdxLocal
        dualSignsLocal
        
        edgeDofsGrouped
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
        function obj = FetiDPSolver(globalMesh, subMeshes, stiffness, forces, tol, dofsNode, boundaryConditions, localToGlobalMaps, useMatrixFree, useEdgeAverage)
            obj.meshCoords          = globalMesh.coord;
            obj.dofsPerNode         = dofsNode;
            obj.localMeshes         = subMeshes;
            obj.localStiffness      = stiffness;
            obj.localForces         = forces;
            obj.nodeTol             = tol;
            obj.numSubdomains       = length(stiffness);
            obj.boundaryConditions  = boundaryConditions;
            obj.useMatrixFree       = useMatrixFree;
            obj.useEdgeAverage      = useEdgeAverage;
            
            obj.primalIdxLocal      = cell(obj.numSubdomains, 1);
            obj.edgeDofsGrouped     = cell(obj.numSubdomains, 1);
            
            obj.extractFetiDofs(localToGlobalMaps);
        end
        
        % -----------------------------------------------------------------
        % 2. MAIN PROBLEM ASSEMBLY
        % -----------------------------------------------------------------
        function [fMatOut, dBar] = assembleProblem(obj)
            allPrimals = unique(vertcat(obj.primalDofsGlobal{:}));
            allDuals   = setdiff(unique(vertcat(obj.dualDofsGlobal{:})), allPrimals);
            
            numPrimals = length(allPrimals);
            numDuals   = length(allDuals);
            
            dBar         = zeros(numDuals, 1);
            sppBase      = sparse(numPrimals, numPrimals);
            BrKrrInvKrp  = sparse(numDuals, numPrimals);
            rhsPrimal    = zeros(numPrimals, 1);
            fDual        = sparse(numDuals, numDuals);
            visitedDuals = zeros(max(allDuals), 1);
            
            if obj.useMatrixFree
                obj.krrFactors = cell(obj.numSubdomains, 1);
                obj.kiiFactors = cell(obj.numSubdomains, 1);
                obj.tdrLocal   = cell(obj.numSubdomains, 1);
            end
            
            for subId = 1:obj.numSubdomains
                [Krr, Krp, Kpr, Kpp, fR, fP, Tdr] = obj.splitLocalMatrices(subId);
                
                if obj.useMatrixFree
                    decKrr = decomposition(Krr);
                    obj.krrFactors{subId} = decKrr;
                    obj.tdrLocal{subId}   = Tdr;
                    
                    rLoc = obj.remDofsLocal{subId};
                    obj.kiiFactors{subId} = decomposition(Krr(1:length(rLoc), 1:length(rLoc)));
                    
                    Urd = decKrr \ Tdr';
                    Urp = decKrr \ Krp;
                else
                    Urd = Krr \ Tdr';
                    Urp = Krr \ Krp;
                end
                
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
                    if ~obj.useMatrixFree
                        fDual(dualRows, dualRows) = fDual(dualRows, dualRows) + dualSigns .* (Urd' * Tdr') .* dualSigns';
                    end
                    dBar(dualRows)               = dBar(dualRows) + dualSigns .* (Urd' * fR);
                    BrKrrInvKrp(dualRows, pRows) = BrKrrInvKrp(dualRows, pRows) + dualSigns .* (Tdr * Urp);
                end
                
                rhsPrimal(pRows) = rhsPrimal(pRows) + (fP - Urp' * fR);
            end
            
            activeDofs    = obj.getActivePrimalDofs(allPrimals);
            obj.sppActive = sppBase(activeDofs, activeDofs);
            obj.brActive  = BrKrrInvKrp(:, activeDofs);
            rhsActive     = rhsPrimal(activeDofs);
            
            dBar = dBar - obj.brActive * (obj.sppActive \ rhsActive);
            
            if obj.useMatrixFree
                fMatOut = @(lambda) obj.applyGlobalF(lambda);
            else
                fMatOut = fDual + obj.brActive * (obj.sppActive \ obj.brActive');
            end
        end
        
        function fOut = applyGlobalF(obj, lambda)
            numDuals = length(lambda);
            fOut     = zeros(numDuals, 1);
            
            for subId = 1:obj.numSubdomains
                dualRows  = obj.dualIdxLocal{subId};
                dualSigns = obj.dualSignsLocal{subId};
                if isempty(dualRows), continue; end
                
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
                if isempty(dualRows), continue; end
                
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
                if isempty(dualRows), continue; end
                
                w    = 1 ./ multiplicity(dualRows);
                rLoc = dualSigns .* r(dualRows);
                
                if obj.useMatrixFree
                    sddLoc = obj.applyLocalSchurVector(subId, w .* rLoc);
                else
                    Sdd = obj.computeLocalSchur(subId);
                    sddLoc = Sdd * (w .* rLoc);
                end
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
                [Krr, Krp, Kpr, Kpp, fR, fP, Tdr] = obj.splitLocalMatrices(subId);
                pRows = obj.getPrimalRows(subId, allPrimals);
                
                if obj.useMatrixFree
                    decKrr = obj.krrFactors{subId};
                    Urp      = decKrr \ Krp;
                    Urd      = decKrr \ Tdr';
                    KrrInvFr = decKrr \ fR;
                else
                    Urp      = Krr \ Krp;
                    Urd      = Krr \ Tdr';
                    KrrInvFr = Krr \ fR;
                end
                
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
            sppAct     = sppBase(activeDofs, activeDofs);
            rhsActive  = rhsPrimal(activeDofs);
            
            uP = zeros(numP, 1);
            uP(activeDofs) = sppAct \ rhsActive;
            
            uFull = zeros(numNodes * obj.dofsPerNode, 1);
            
            for subId = 1:obj.numSubdomains
                pGlobal = obj.primalDofsGlobal{subId};
                dGlobal = obj.dualDofsGlobal{subId};
                rGlobal = obj.remDofsGlobal{subId};
                
                pLocal = obj.primalDofsLocal{subId};
                dLocal = obj.dualDofsLocal{subId};
                rLocal = obj.remDofsLocal{subId};
                pRows  = obj.getPrimalRows(subId, allPrimals);
                
                dualRows  = obj.dualIdxLocal{subId};
                dualSigns = obj.dualSignsLocal{subId};
                lambdaLoc = dualSigns .* lambdaSol(dualRows);
                
                uPLoc   = uP(pRows);
                uRemLoc = krrInvFrCell{subId} - urpCell{subId} * uPLoc - urdCell{subId} * lambdaLoc;
                
                numLocalDofs = length(pLocal) + length(dLocal) + length(rLocal);
                uLocalTilde = zeros(numLocalDofs, 1);
                uLocalTilde(pLocal) = uPLoc;
                uLocalTilde([rLocal; dLocal]) = uRemLoc;
                
                if obj.useEdgeAverage
                    edgeGroups = obj.edgeDofsGrouped{subId};
                    uLocalPhys = obj.applyEdgeAverageForward(uLocalTilde, edgeGroups);
                else
                    uLocalPhys = uLocalTilde;
                end
                
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
            
            pCoords = obj.meshCoords(primalNodes, :);
            dCoords = obj.meshCoords(dualNodes, :);
            rCoords = obj.meshCoords(remNodes, :);
            
            figure('Name', 'Distribución de Nodos FETI-DP', 'Color', 'w');
            hold on; axis equal;
            for i = 1:obj.numSubdomains
                patch('Faces', obj.localMeshes{i}.connec, 'Vertices', obj.localMeshes{i}.coord, ...
                    'FaceColor', 'none', 'EdgeColor', [0.8 0.8 0.8], 'LineWidth', 1.5, 'HandleVisibility', 'off');
            end
            scatter(rCoords(:, 1), rCoords(:, 2), 20, [0.5 0.5 0.5], 'filled', 'DisplayName', 'Interior');
            scatter(dCoords(:, 1), dCoords(:, 2), 40, 'b', 'filled', 'DisplayName', 'Interface (Dual)');
            scatter(pCoords(:, 1), pCoords(:, 2), 80, 'r', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1, 'DisplayName', 'Primal (Esquinas/Bordes)');
            legend('Location', 'bestoutside'); grid on; hold off;
        end
    end
    
    methods (Access = private)
        
        % -----------------------------------------------------------------
        % 6. TOPOLOGICAL SETUP
        % -----------------------------------------------------------------
        function extractFetiDofs(obj, localToGlobalMaps)
            nSub = obj.numSubdomains;
            
            % Preallocate cell arrays for DOFs
            pGlobal = cell(nSub, 1); 
            dGlobal = cell(nSub, 1); 
            rGlobal = cell(nSub, 1);
            pLocal  = cell(nSub, 1); 
            dLocal  = cell(nSub, 1); 
            rLocal  = cell(nSub, 1);
            
            % 1. GLOBAL INTERFACE MAPPING & SUBDOMAIN MASK
            numGlobalNodes   = size(obj.meshCoords, 1);
            nodeMultiplicity = zeros(numGlobalNodes, 1);
            nodeSubMask      = sparse(numGlobalNodes, nSub);
            
            for i = 1:nSub
                gN = localToGlobalMaps{i};
                nodeMultiplicity(gN) = nodeMultiplicity(gN) + 1;
                nodeSubMask(gN, i) = 1; 
            end
            
            % 2. GEOMETRIC CLASSIFICATION PER SUBDOMAIN
            for i = 1:nSub
                gNodes = localToGlobalMaps{i};
                nNodes = length(gNodes);
                
                % --- A. Initial Primal Nodes (Corners) ---
                boundaryMeshes = obj.localMeshes{i}.createBoundaryMesh();
                numBoundaryMeshes = length(boundaryMeshes);
                primalNodesPerMesh = cell(numBoundaryMeshes, 1);
                
                for b = 1:numBoundaryMeshes
                    boundaryConnectivity = boundaryMeshes{b}.globalConnec(:);                    
                    nodeOccurrences = accumarray(boundaryConnectivity, 1);
                    primalNodesPerMesh{b} = find(nodeOccurrences == 1);
                end                
                idxPrimal = unique(vertcat(primalNodesPerMesh{:}));
                
                % --- B. Dual Nodes and Edge Grouping ---
                if obj.useEdgeAverage
                    isShared = nodeMultiplicity(gNodes) > 1;
                    idxDualOrig = setdiff(find(isShared), idxPrimal);
                    gDual = gNodes(idxDualOrig);
                    
                    % Identify unique edges shared by the exact same set of subdomains
                    [~, ~, edgeIdx] = unique(full(nodeSubMask(gDual, :)), 'rows');
                    numEdges = max(edgeIdx);
                    
                    edgeDofsLocalSub = cell(numEdges, 1);
                    idxDualFinal = [];
                    
                    for e = 1:numEdges
                        localNodesInEdge = idxDualOrig(edgeIdx == e);
                        [~, sortIdx] = sort(gNodes(localNodesInEdge));
                        localNodesInEdge = localNodesInEdge(sortIdx);
                        
                        if length(localNodesInEdge) > 1
                            edgeDofsLocalSub{e} = obj.nodesToDofs(localNodesInEdge);
                            
                            % FETI-DP constraint: First edge node to primal, rest to dual
                            idxPrimal    = [idxPrimal; localNodesInEdge(1)];
                            idxDualFinal = [idxDualFinal; localNodesInEdge(2:end)];
                        else
                            % Isolated interface node becomes standard dual
                            idxDualFinal = [idxDualFinal; localNodesInEdge(1)];
                        end
                    end
                    
                    idxDualFinal = unique(idxDualFinal);
                    
                    % Store active edge DOFs filtering empty cells
                    obj.edgeDofsGrouped{i} = edgeDofsLocalSub(~cellfun('isempty', edgeDofsLocalSub));
                    
                else
                    % --- B. Standard Dual Nodes ---
                    isShared = nodeMultiplicity(gNodes) > 1;
                    idxDualFinal = setdiff(find(isShared), idxPrimal);
                end
                
                % --- C. Remaining Nodes (Interior) ---
                idxRem = setdiff((1:nNodes)', [idxPrimal; idxDualFinal]);
                
                % 3. NODE TO DOF CONVERSION
                pGlobal{i} = obj.nodesToDofs(gNodes(idxPrimal));
                dGlobal{i} = obj.nodesToDofs(gNodes(idxDualFinal));
                rGlobal{i} = obj.nodesToDofs(gNodes(idxRem));
                
                pLocal{i} = obj.nodesToDofs(idxPrimal);
                dLocal{i} = obj.nodesToDofs(idxDualFinal);
                rLocal{i} = obj.nodesToDofs(idxRem);
            end
            
            % 4. CLASS PROPERTY ASSIGNMENT
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
            
            if obj.useEdgeAverage
                edgeGroups = obj.edgeDofsGrouped{subId};
                [kMat, fVec] = obj.applyEdgeAverageTransformation(kMat, fVec, edgeGroups);
            end
            
            remLoc = [rLoc; dLoc];
            Krr = kMat(remLoc, remLoc);
            Krp = kMat(remLoc, pLoc);
            Kpr = kMat(pLoc, remLoc);
            Kpp = kMat(pLoc, pLoc);
            
            fR  = fVec(remLoc);
            fP  = fVec(pLoc);
            
            numD = length(dLoc); numR = length(remLoc);
            Tdr  = sparse((1:numD)', (length(rLoc) + 1 : numR)', ones(numD, 1), numD, numR);
        end
        
        function T = buildEdgeTransform(obj, numDofs, edgeGroups)
            T = speye(numDofs);
            for e = 1:length(edgeGroups)
                eDofs = edgeGroups{e};
                for d = 1:obj.dofsPerNode
                    dimDofs = eDofs(d:obj.dofsPerNode:end);
                    T(dimDofs, dimDofs(1)) = 1;
                    T(dimDofs(1), dimDofs(2:end)) = -1;
                end
            end
        end
        
        function [kMod, fMod] = applyEdgeAverageTransformation(obj, K, f, edgeGroups)
            T = obj.buildEdgeTransform(size(K, 1), edgeGroups);
            kMod = T' * K * T; fMod = T' * f;
        end
        
        function uPhys = applyEdgeAverageForward(obj, uTransformed, edgeGroups)
            T = obj.buildEdgeTransform(length(uTransformed), edgeGroups);
            uPhys = T * uTransformed;
        end
        
        % -----------------------------------------------------------------
        % 8. MATHEMATICAL OPERATIONS
        % -----------------------------------------------------------------
        function sddX = applyLocalSchurVector(obj, subId, x)
            kMat = obj.localStiffness{subId};
            rLoc = obj.remDofsLocal{subId};
            dLoc = obj.dualDofsLocal{subId};
            if obj.useEdgeAverage
                edgeGroups = obj.edgeDofsGrouped{subId};
                [kMat, ~]  = obj.applyEdgeAverageTransformation(kMat, zeros(size(kMat,1),1), edgeGroups);
            end
            Kid = kMat(rLoc, dLoc); Kdi = kMat(dLoc, rLoc); Kdd = kMat(dLoc, dLoc);
            sddX = Kdd * x - Kdi * (obj.kiiFactors{subId} \ (Kid * x));
        end
        
        function Sdd = computeLocalSchur(obj, subId)
            kMat = obj.localStiffness{subId};
            rLoc = obj.remDofsLocal{subId};
            dLoc = obj.dualDofsLocal{subId};
            if obj.useEdgeAverage
                edgeGroups = obj.edgeDofsGrouped{subId};
                [kMat, ~]  = obj.applyEdgeAverageTransformation(kMat, zeros(size(kMat,1),1), edgeGroups);
            end
            Sdd = kMat(dLoc, dLoc) - kMat(dLoc, rLoc) * (kMat(rLoc, rLoc) \ kMat(rLoc, dLoc));
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