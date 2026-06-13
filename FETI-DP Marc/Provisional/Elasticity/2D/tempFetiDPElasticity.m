classdef tempFetiDPElasticity < handle
    % FETI-DP algebraic solver for 2D linear elasticity on a domain decomposition.
    
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
    end
    
    methods (Access = public)
        
        % -----------------------------------------------------------------
        % 1. CONSTRUCTOR
        % -----------------------------------------------------------------
        function obj = tempFetiDPElasticity(globalMesh, subMeshes, stiffness, forces, tol, dofsNode, boundaryConditions, localToGlobalMaps)
            obj.meshCoords          = globalMesh.coord;
            obj.dofsPerNode         = dofsNode;
            obj.localMeshes         = subMeshes;
            obj.localStiffness      = stiffness;
            obj.localForces         = forces;
            obj.nodeTol             = tol;
            obj.numSubdomains       = length(stiffness);
            obj.boundaryConditions  = boundaryConditions;
            obj.primalIdxLocal      = cell(obj.numSubdomains, 1);
            
            obj.extractFetiDofs(localToGlobalMaps);
        end
        
        % -----------------------------------------------------------------
        % 2. MAIN PROBLEM ASSEMBLY
        % -----------------------------------------------------------------
        function [fMat, dBar] = assembleProblem(obj)
            allPrimals = unique(vertcat(obj.primalDofsGlobal{:}));
            allDuals   = setdiff(unique(vertcat(obj.dualDofsGlobal{:})), allPrimals);
            
            numPrimals = length(allPrimals);
            numDuals   = length(allDuals);
            
            fDual       = sparse(numDuals, numDuals);
            dBar        = zeros(numDuals, 1);
            SPP         = sparse(numPrimals, numPrimals);
            BrKrrInvKrp = sparse(numDuals, numPrimals);
            rhsPrimal   = zeros(numPrimals, 1);
            
            visitedDuals = zeros(max(allDuals), 1);
            
            for subId = 1:obj.numSubdomains
                [Krr, Krp, Kpr, Kpp, fR, fP, Tdr] = obj.splitLocalMatrices(subId);
                
                Urd = Krr \ Tdr';
                Urp = Krr \ Krp;
                
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
                SPP(pRows, pRows) = SPP(pRows, pRows) + sppLoc;
                
                if ~isempty(dGlobal)
                    fDual(dualRows, dualRows)       = fDual(dualRows, dualRows) + dualSigns .* (Urd' * Tdr') .* dualSigns';
                    dBar(dualRows)                  = dBar(dualRows) + dualSigns .* (Urd' * fR);
                    BrKrrInvKrp(dualRows, pRows)    = BrKrrInvKrp(dualRows, pRows) + dualSigns .* (Tdr * Urp);
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
                
                Wsigned = diag(w .* dualSigns);
                M(dualRows, dualRows) = M(dualRows, dualRows) + Wsigned * Sdd * Wsigned;
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
                
                Sdd  = obj.computeLocalSchur(subId);
                w    = 1 ./ multiplicity(dualRows);
                rLoc = dualSigns .* r(dualRows);
                
                z(dualRows) = z(dualRows) + dualSigns .* (w .* (Sdd * (w .* rLoc)));
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
                
                dualRows  = obj.dualIdxLocal{subId};
                dualSigns = obj.dualSignsLocal{subId};
                lambdaLoc = dualSigns .* lambdaSol(dualRows);
                term      = Tdr' * lambdaLoc - fR;
                rhsPrimal(pRows) = rhsPrimal(pRows) + (fP + Urp' * term);
            end
            
            activeDofs = obj.getActivePrimalDofs(allPrimals);
            sppActive  = SPP(activeDofs, activeDofs);
            rhsActive  = rhsPrimal(activeDofs);
            
            uP = zeros(numP, 1);
            uP(activeDofs) = sppActive \ rhsActive;
            
            uFull = zeros(numNodes * obj.dofsPerNode, 1);
            
            for subId = 1:obj.numSubdomains
                pGlobal = obj.primalDofsGlobal{subId};
                dGlobal = obj.dualDofsGlobal{subId};
                rGlobal = obj.remDofsGlobal{subId};
                pRows   = obj.getPrimalRows(subId, allPrimals);
                
                dualRows  = obj.dualIdxLocal{subId};
                dualSigns = obj.dualSignsLocal{subId};
                lambdaLoc = dualSigns .* lambdaSol(dualRows);
                
                uPLoc   = uP(pRows);
                uRemLoc = KrrInvFrCell{subId} - UrpCell{subId} * uPLoc - UrdCell{subId} * lambdaLoc;
                
                uFull(pGlobal)            = uPLoc;
                uFull([rGlobal; dGlobal]) = uRemLoc;
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
            scatter(pCoords(:, 1), pCoords(:, 2), 80, 'r', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1, 'DisplayName', 'Corners (Primal)');
            
            legend('Location', 'bestoutside');
            title('FETI-DP Node Distribution');
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
            tol  = obj.nodeTol;
            
            % Preallocate cell arrays for DOFs
            pGlobal = cell(nSub, 1); 
            dGlobal = cell(nSub, 1); 
            rGlobal = cell(nSub, 1);
            pLocal  = cell(nSub, 1); 
            dLocal  = cell(nSub, 1); 
            rLocal  = cell(nSub, 1);
            
            % 1. GLOBAL INTERFACE MAPPING
            numGlobalNodes   = size(obj.meshCoords, 1);
            nodeMultiplicity = zeros(numGlobalNodes, 1);
            
            for i = 1:nSub
                gN = localToGlobalMaps{i};
                nodeMultiplicity(gN) = nodeMultiplicity(gN) + 1;
            end
            
            % 2. GEOMETRIC CLASSIFICATION PER SUBDOMAIN
            for i = 1:nSub
                gNodes = localToGlobalMaps{i};
                nNodes = length(gNodes);
                
                % --- A. Primal Nodes ---
                boundaryMeshes = obj.localMeshes{i}.createBoundaryMesh();
                numBoundaryMeshes = length(boundaryMeshes);
                primalNodesPerMesh = cell(numBoundaryMeshes, 1);
                
                for b = 1:numBoundaryMeshes
                    boundaryConnectivity = boundaryMeshes{b}.globalConnec(:);                    
                    nodeOccurrences = accumarray(boundaryConnectivity, 1);
                    primalNodesPerMesh{b} = find(nodeOccurrences == 1);
                end                
                idxPrimal = unique(vertcat(primalNodesPerMesh{:}));
                
                % --- Dual Nodes ---
                isShared = nodeMultiplicity(gNodes) > 1;
                idxDual = setdiff(find(isShared), idxPrimal);
                
                % --- Remaining Nodes (Interior) ---
                idxRem = setdiff((1:nNodes)', [idxPrimal; idxDual]);
                
                % 3. NODE TO DOF CONVERSION
                pGlobal{i} = obj.nodesToDofs(gNodes(idxPrimal));
                dGlobal{i} = obj.nodesToDofs(gNodes(idxDual));
                rGlobal{i} = obj.nodesToDofs(gNodes(idxRem));
                
                pLocal{i} = obj.nodesToDofs(idxPrimal);
                dLocal{i} = obj.nodesToDofs(idxDual);
                rLocal{i} = obj.nodesToDofs(idxRem);
            end
            
            obj.primalDofsGlobal = pGlobal; obj.primalDofsLocal = pLocal;
            obj.dualDofsGlobal   = dGlobal; obj.dualDofsLocal   = dLocal;
            obj.remDofsGlobal    = rGlobal; obj.remDofsLocal    = rLocal;
        end

        function dofs = nodesToDofs(obj, nodes)
            nodes = nodes(:);
            base = (nodes - 1) * obj.dofsPerNode;
            dofs = base + (1:obj.dofsPerNode);
            dofs = dofs.'; 
            dofs = dofs(:); 
        end
        
        % -----------------------------------------------------------------
        % 7. ALGEBRAIC TRANSFORMATIONS
        % -----------------------------------------------------------------
        function [Krr, Krp, Kpr, Kpp, fR, fP, Tdr] = splitLocalMatrices(obj, subId)
            kMat   = obj.localStiffness{subId};
            fVec   = obj.localForces{subId};
            
            pLoc   = obj.primalDofsLocal{subId};
            dLoc   = obj.dualDofsLocal{subId};
            rLoc   = obj.remDofsLocal{subId};
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
        
        function activeIdx = getActivePrimalDofs(~, allPrimals)
            activeIdx = (1:length(allPrimals));
        end
        
        
    end
end
