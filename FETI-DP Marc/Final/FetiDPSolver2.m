classdef FetiDPSolver2 < handle
    % -----------------------------------------------------------------
    % FETI-DP ALGEBRAIC SOLVER (INDEX-BASED CHAIN TOPOLOGY & BOUNDARY FIX)
    % -----------------------------------------------------------------
    
    properties (Access = public)
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
        
        lambdaGlobalIdx
        lambdaLocalIdx
        lambdaSigns
        dualWeightLocal
        numGlobalLambdas
        
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
        function obj = FetiDPSolver2(globalMesh, subMeshes, stiffness, forces, tol, dofsNode, boundaryConditions, localToGlobalMaps, useMatrixFree, useEdgeAverage)
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
            numPrimals = length(allPrimals);
            numDuals   = obj.numGlobalLambdas;
            
            dBar         = zeros(numDuals, 1);
            sppBase      = sparse(numPrimals, numPrimals);
            BrKrrInvKrp  = sparse(numDuals, numPrimals);
            rhsPrimal    = zeros(numPrimals, 1);
            fDual        = sparse(numDuals, numDuals);
            
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
                
                sppLoc = Kpp - Kpr * Urp;
                sppBase(pRows, pRows) = sppBase(pRows, pRows) + sppLoc;
                
                gIdx = obj.lambdaGlobalIdx{subId};
                lIdx = obj.lambdaLocalIdx{subId};
                sgn  = obj.lambdaSigns{subId};
                
                loc_dBar = Urd' * fR;
                dBar(gIdx) = dBar(gIdx) + sgn .* loc_dBar(lIdx);
                
                loc_Br = Tdr * Urp;
                BrKrrInvKrp(gIdx, pRows) = BrKrrInvKrp(gIdx, pRows) + sgn .* loc_Br(lIdx, :);
                
                if ~obj.useMatrixFree
                    loc_F = Urd' * Tdr';
                    for k1 = 1:length(gIdx)
                        for k2 = 1:length(gIdx)
                            fDual(gIdx(k1), gIdx(k2)) = fDual(gIdx(k1), gIdx(k2)) + sgn(k1) * loc_F(lIdx(k1), lIdx(k2)) * sgn(k2);
                        end
                    end
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
                gIdx = obj.lambdaGlobalIdx{subId};
                lIdx = obj.lambdaLocalIdx{subId};
                sgn  = obj.lambdaSigns{subId};
                                
                numLocalDuals = length(obj.dualDofsLocal{subId});
                lambdaLoc     = accumarray(lIdx, sgn .* lambda(gIdx), [numLocalDuals, 1]);
                
                Tdr     = obj.tdrLocal{subId};
                decKrr  = obj.krrFactors{subId};
                termLoc = Tdr * (decKrr \ (Tdr' * lambdaLoc));
                
                fOut(gIdx) = fOut(gIdx) + sgn .* termLoc(lIdx);
            end
            
            primalCorrection = obj.brActive * (obj.sppActive \ (obj.brActive' * lambda));
            fOut = fOut + primalCorrection;
        end
        
        % -----------------------------------------------------------------
        % 3. SOLVER TOOLS (PRECONDITIONER)
        % -----------------------------------------------------------------
        function M = buildPrecondMatrix(obj)
            numDuals = obj.numGlobalLambdas;
            M = zeros(numDuals);
            
            for subId = 1:obj.numSubdomains
                gIdx = obj.lambdaGlobalIdx{subId};
                lIdx = obj.lambdaLocalIdx{subId};
                sgn  = obj.lambdaSigns{subId};
                W    = obj.dualWeightLocal{subId};
                                
                Sdd = obj.computeLocalSchur(subId);
                
                for k1 = 1:length(gIdx)
                    for k2 = 1:length(gIdx)
                        M(gIdx(k1), gIdx(k2)) = M(gIdx(k1), gIdx(k2)) + ...
                            sgn(k1) * W(lIdx(k1)) * Sdd(lIdx(k1), lIdx(k2)) * W(lIdx(k2)) * sgn(k2);
                    end
                end
            end
        end
        
        function z = applyDirichletPrecond(obj, r)
            numDuals = obj.numGlobalLambdas;
            z = zeros(numDuals, 1);
            
            for subId = 1:obj.numSubdomains
                gIdx = obj.lambdaGlobalIdx{subId};
                lIdx = obj.lambdaLocalIdx{subId};
                sgn  = obj.lambdaSigns{subId};
                W    = obj.dualWeightLocal{subId};
                                
                numLocalDuals = length(obj.dualDofsLocal{subId});
                rLoc = accumarray(lIdx, sgn .* r(gIdx), [numLocalDuals, 1]);
                
                if obj.useMatrixFree
                    sddLoc = obj.applyLocalSchurVector(subId, W .* rLoc);
                else
                    Sdd = obj.computeLocalSchur(subId);
                    sddLoc = Sdd * (W .* rLoc);
                end
                
                z(gIdx) = z(gIdx) + sgn .* W(lIdx) .* sddLoc(lIdx);
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
                    decKrr   = obj.krrFactors{subId};
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
                
                gIdx = obj.lambdaGlobalIdx{subId};
                lIdx = obj.lambdaLocalIdx{subId};
                sgn  = obj.lambdaSigns{subId};
                numLocalDuals = length(obj.dualDofsLocal{subId});
                
                if ~isempty(gIdx)
                    lambdaLoc = accumarray(lIdx, sgn .* lambdaSol(gIdx), [numLocalDuals, 1]);
                else
                    lambdaLoc = zeros(numLocalDuals, 1);
                end
                
                term = Tdr' * lambdaLoc - fR;
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
                
                gIdx = obj.lambdaGlobalIdx{subId};
                lIdx = obj.lambdaLocalIdx{subId};
                sgn  = obj.lambdaSigns{subId};
                numLocalDuals = length(dLocal);
                
                if ~isempty(gIdx)
                    lambdaLoc = accumarray(lIdx, sgn .* lambdaSol(gIdx), [numLocalDuals, 1]);
                else
                    lambdaLoc = zeros(numLocalDuals, 1);
                end
                
                uPLoc   = uP(pRows);
                uRemLoc = krrInvFrCell{subId} - urpCell{subId} * uPLoc - urdCell{subId} * lambdaLoc;
                
                numLocalDofs = length(pLocal) + length(dLocal) + length(rLocal);
                uLocalTilde  = zeros(numLocalDofs, 1);
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
            
            dim = size(obj.meshCoords, 2);
            
            figure('Name', 'Nodes Distribution FETI-DP', 'Color', 'w');
            hold on; axis equal;
            
            for i = 1:obj.numSubdomains
                patch('Faces', obj.localMeshes{i}.connec, 'Vertices', obj.localMeshes{i}.coord, ...
                    'FaceColor', 'none', 'EdgeColor', [0.8 0.8 0.8], 'LineWidth', 1.5, 'HandleVisibility', 'off');
            end
            
            if dim == 2
                scatter(rCoords(:, 1), rCoords(:, 2), 20, [0.5 0.5 0.5], 'filled', 'DisplayName', 'Remaining/Interior');
                scatter(dCoords(:, 1), dCoords(:, 2), 40, 'b', 'filled', 'DisplayName', 'Interface (Dual)');
                scatter(pCoords(:, 1), pCoords(:, 2), 80, 'r', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1, 'DisplayName', 'Corner (Primal)');
            elseif dim == 3
                scatter3(rCoords(:, 1), rCoords(:, 2), rCoords(:, 3), 20, [0.5 0.5 0.5], 'filled', 'DisplayName', 'Interior');
                scatter3(dCoords(:, 1), dCoords(:, 2), dCoords(:, 3), 40, 'b', 'filled', 'DisplayName', 'Interface (Dual)');
                scatter3(pCoords(:, 1), pCoords(:, 2), pCoords(:, 3), 80, 'r', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1, 'DisplayName', 'Primal');
                zlabel('Z');
                view(3);
            end
            axis off;
            legend('Location', 'east', 'FontSize', 14);
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
            dim              = size(obj.meshCoords, 2); 
            nodeMultiplicity = zeros(numGlobalNodes, 1);
            nodeSubMask      = sparse(numGlobalNodes, nSub);
            
            for i = 1:nSub
                gN = localToGlobalMaps{i};
                nodeMultiplicity(gN) = nodeMultiplicity(gN) + 1;
                nodeSubMask(gN, i) = 1;  
            end
            
            % --- DETECT GLOBAL BOUNDARY (only used for primal corner detection) ---
            tolCoords = 1e-6;
            minC = min(obj.meshCoords);
            maxC = max(obj.meshCoords);
            
            if dim == 2
                isOnGlobalBoundary = (abs(obj.meshCoords(:,1) - minC(1)) < tolCoords) | ...
                                     (abs(obj.meshCoords(:,1) - maxC(1)) < tolCoords) | ...
                                     (abs(obj.meshCoords(:,2) - minC(2)) < tolCoords) | ...
                                     (abs(obj.meshCoords(:,2) - maxC(2)) < tolCoords);
            elseif dim == 3
                isOnGlobalBoundary = (abs(obj.meshCoords(:,1) - minC(1)) < tolCoords) | ...
                                     (abs(obj.meshCoords(:,1) - maxC(1)) < tolCoords) | ...
                                     (abs(obj.meshCoords(:,2) - minC(2)) < tolCoords) | ...
                                     (abs(obj.meshCoords(:,2) - maxC(2)) < tolCoords) | ...
                                     (abs(obj.meshCoords(:,3) - minC(3)) < tolCoords) | ...
                                     (abs(obj.meshCoords(:,3) - maxC(3)) < tolCoords);
            end
            
            for i = 1:nSub
                gNodes = localToGlobalMaps{i};
                nNodes = length(gNodes);
                
                boundaryMeshes    = obj.localMeshes{i}.createBoundaryMesh();
                numBoundaryMeshes = length(boundaryMeshes);
                
                if dim == 2
                    % 2D: nodo primal = aparece solo 1 vez en el contorno (esquina)
                    primalNodesPerMesh = cell(numBoundaryMeshes, 1);
                    for b = 1:numBoundaryMeshes
                        boundaryConnectivity   = boundaryMeshes{b}.globalConnec(:);
                        nodeOccurrences        = accumarray(boundaryConnectivity, 1);
                        primalNodesPerMesh{b}  = find(nodeOccurrences == 1);
                    end
                    idxPrimal = unique(vertcat(primalNodesPerMesh{:}));
                else
                    % 3D: nodo primal = aparece en >= 3 caras del boundary local
                    % Acumular TODOS los nodos de TODAS las caras de una vez.
                    allBoundaryNodes = [];
                    for bb = 1:numBoundaryMeshes
                        nodesInFace      = unique(boundaryMeshes{bb}.globalConnec(:));
                        allBoundaryNodes = [allBoundaryNodes; nodesInFace];
                    end
                    nodeOccurrences = accumarray(allBoundaryNodes, 1);
                    idxPrimal = find(nodeOccurrences >= 3);
                    idxPrimal = idxPrimal(:);
                end
                
                if obj.useEdgeAverage
                    isShared = nodeMultiplicity(gNodes) > 1;
                    idxDualOrig = setdiff(find(isShared), idxPrimal);
                    gDual = gNodes(idxDualOrig);
                    
                    % --- AGRUPAR POR CONJUNTO DE SUBDOMINIOS COMPARTIDOS ---
                    % Cada "arista" es el conjunto de nodos que comparten exactamente
                    % los mismos subdominios. La clave es la fila de nodeSubMask,
                    % sin usar coordenadas.
                    [~, ~, edgeIdx] = unique(full(nodeSubMask(gDual, :)), 'rows');
                    numEdges = max(edgeIdx);
                    
                    edgeDofsLocalSub = cell(numEdges, 1);
                    idxDualFinal = [];
                    
                    for e = 1:numEdges
                        localNodesInEdge = idxDualOrig(edgeIdx == e);
                        % Ordenar por índice global para determinismo
                        [~, sortIdx] = sort(gNodes(localNodesInEdge));
                        localNodesInEdge = localNodesInEdge(sortIdx);
                        
                        % Siempre asignar un primal al primer nodo del grupo,
                        % incluso si la arista tiene un solo nodo.
                        edgeDofsLocalSub{e} = obj.nodesToDofs(localNodesInEdge);
                        idxPrimal    = [idxPrimal; localNodesInEdge(1)];
                        if length(localNodesInEdge) > 1
                            idxDualFinal = [idxDualFinal; localNodesInEdge(2:end)];
                        end
                        % Si solo hay un nodo en la arista, queda como primal
                        % (no se añade a idxDualFinal).
                    end
                    idxDualFinal = unique(idxDualFinal);
                    obj.edgeDofsGrouped{i} = edgeDofsLocalSub(~cellfun('isempty', edgeDofsLocalSub));
                    
                else
                    isShared = nodeMultiplicity(gNodes) > 1;
                    idxDualFinal = setdiff(find(isShared), idxPrimal);
                end
                
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
            
            % -------------------------------------------------------------
            % 6.1. CHAIN TOPOLOGY (INDEX GENERATION)
            % -------------------------------------------------------------
            allDualMeshDofs  = unique(vertcat(obj.dualDofsGlobal{:}));
            numPhysicalDuals = length(allDualMeshDofs);
            
            if numPhysicalDuals == 0
                obj.numGlobalLambdas = 0;
                return;
            end
            
            maxMeshDof = max(allDualMeshDofs);
            meshDofToIdx = zeros(maxMeshDof, 1);
            meshDofToIdx(allDualMeshDofs) = 1:numPhysicalDuals;
            
            subdomainsSharing = cell(numPhysicalDuals, 1);
            localIdxSharing   = cell(numPhysicalDuals, 1);
            
            for subId = 1:nSub
                dGlob = obj.dualDofsGlobal{subId};
                if isempty(dGlob), continue; end
                
                idx = meshDofToIdx(dGlob);
                for k = 1:length(idx)
                    physId = idx(k);
                    subdomainsSharing{physId}(end+1) = subId;
                    localIdxSharing{physId}(end+1)   = k;
                end
            end
            
            obj.lambdaGlobalIdx = cell(nSub, 1);
            obj.lambdaLocalIdx  = cell(nSub, 1);
            obj.lambdaSigns     = cell(nSub, 1);
            obj.dualWeightLocal = cell(nSub, 1);
            
            for subId = 1:nSub
                numLocalDuals = length(obj.dualDofsLocal{subId});
                obj.dualWeightLocal{subId} = ones(numLocalDuals, 1);
            end
            
            lambdaCounter = 1;
            
            for physId = 1:numPhysicalDuals
                subs = subdomainsSharing{physId};
                locs = localIdxSharing{physId};
                M    = length(subs);
                
                for i = 1:M
                    obj.dualWeightLocal{subs(i)}(locs(i)) = 1.0 / M;
                end
                
                for i = 1:M-1
                    subA = subs(i);   locA = locs(i);
                    subB = subs(i+1); locB = locs(i+1);
                    
                    obj.lambdaGlobalIdx{subA}(end+1, 1) = lambdaCounter;
                    obj.lambdaLocalIdx{subA}(end+1, 1)  = locA;
                    obj.lambdaSigns{subA}(end+1, 1)     = 1;
                    
                    obj.lambdaGlobalIdx{subB}(end+1, 1) = lambdaCounter;
                    obj.lambdaLocalIdx{subB}(end+1, 1)  = locB;
                    obj.lambdaSigns{subB}(end+1, 1)     = -1;
                    
                    lambdaCounter = lambdaCounter + 1;
                end
            end
            
            obj.numGlobalLambdas = lambdaCounter - 1;
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
            kMod = T' * K * T; 
            fMod = T' * f;
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