classdef FetiDPElasticityFV < handle
    % FETIDPELASTICITY Class for the algebraic logic of the FETI-DP method in 2D Elasticity.
    
    properties (Access = private)
        % Subdomain & Mesh Data
        localStiffness
        localForces
        localMeshes
        meshCoords
        dofsPerNode
        nodeTol
        numSubdomains
        dirichletDofs
        
        % Global DoF Mappings
        primalDofsGlobal
        dualDofsGlobal
        remDofsGlobal
        
        % Local DoF Mappings
        primalDofsLocal
        dualDofsLocal
        remDofsLocal
        edgeNodesLocal
        Tlocal
        
        % Preconditioner & Interface Assembly Data
        dualIdxLocal    
        dualSignsLocal  
        primalIdxLocal
        localSchurBlocks % Precomputed Schur complements for the preconditioner
    end
    
    % =========================================================
    % PUBLIC METHODS: INITIALIZATION & MAIN SOLVER
    % =========================================================
    methods (Access = public)
        
        function obj = FetiDPElasticityFV(globalMesh, subMeshes, stiffness, forces, tol, dofsNode, dirDofs)
            obj.meshCoords     = globalMesh.coord;
            obj.dofsPerNode    = dofsNode;
            obj.localMeshes    = subMeshes;
            obj.localStiffness = stiffness;
            obj.localForces    = forces;
            obj.nodeTol        = tol;
            obj.numSubdomains  = length(stiffness);
            obj.dirichletDofs  = dirDofs;
            
            obj.extractFetiDofs();
            obj.applyEdgeTransformation();
        end
        
        function [fMat, dBar] = assembleProblem(obj)
            allPrimals = unique(vertcat(obj.primalDofsGlobal{:}));
            allDuals   = setdiff(unique(vertcat(obj.dualDofsGlobal{:})), allPrimals);
            
            numPrimals = length(allPrimals);
            numDuals   = length(allDuals);
            
            fDual        = sparse(numDuals, numDuals);
            dBar         = zeros(numDuals, 1);
            SPP          = sparse(numPrimals, numPrimals);
            BrKrrInvKrp  = sparse(numDuals, numPrimals);
            rhsPrimal    = zeros(numPrimals, 1);
            
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
                
                isFirst             = ~visitedDuals(dGlobal);
                dualSigns           = ones(length(dGlobal), 1);
                dualSigns(~isFirst) = -1;
                visitedDuals(dGlobal) = true;
                
                obj.dualIdxLocal{subId}   = dualRows;
                obj.dualSignsLocal{subId} = dualSigns;
                
                SppLoc            = Kpp - Kpr * Urp;
                SPP(pRows, pRows) = SPP(pRows, pRows) + SppLoc;
                
                if ~isempty(dGlobal)
                    fDual(dualRows, dualRows)    = fDual(dualRows, dualRows) + dualSigns .* (Urd' * Tdr') .* dualSigns';
                    dBar(dualRows)               = dBar(dualRows) + dualSigns .* (Urd' * fR);
                    BrKrrInvKrp(dualRows, pRows) = BrKrrInvKrp(dualRows, pRows) + dualSigns .* (Tdr * Urp);
                end
                
                rhsPrimal(pRows) = rhsPrimal(pRows) + (fP - Urp' * fR);
            end
            
            activeDofs = obj.getActivePrimalDofs(allPrimals);
            SPPActive  = SPP(activeDofs, activeDofs);
            BrActive   = BrKrrInvKrp(:, activeDofs);
            rhsActive  = rhsPrimal(activeDofs);
            
            fMat = fDual + BrActive * (SPPActive \ BrActive');
            dBar = dBar  - BrActive * (SPPActive \ rhsActive);
            
            % Precompute local Schur blocks for fast preconditioner application
            obj.localSchurBlocks = cell(obj.numSubdomains, 1);
            for subId = 1:obj.numSubdomains
                if ~isempty(obj.dualIdxLocal{subId})
                    obj.localSchurBlocks{subId} = obj.computeLocalSchur(subId);
                end
            end
        end
        
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
                
                pRows = obj.primalIdxLocal{subId}; 
                
                Urp      = Krr \ Krp;
                Urd      = Krr \ Tdr';
                KrrInvFr = Krr \ fR;
                
                UrpCell{subId}      = Urp;
                UrdCell{subId}      = Urd;
                KrrInvFrCell{subId} = KrrInvFr;
                
                sppLoc            = Kpp - Kpr * Urp;
                SPP(pRows, pRows) = SPP(pRows, pRows) + sppLoc;
                
                dualRows  = obj.dualIdxLocal{subId};
                dualSigns = obj.dualSignsLocal{subId};
                
                lambdaLoc        = dualSigns .* lambdaSol(dualRows);
                term             = Tdr' * lambdaLoc - fR;
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
                
                pRows     = obj.primalIdxLocal{subId};
                dualRows  = obj.dualIdxLocal{subId};
                dualSigns = obj.dualSignsLocal{subId};
                lambdaLoc = dualSigns .* lambdaSol(dualRows);
                
                uPLoc   = uP(pRows); 
                uRemLoc = KrrInvFrCell{subId} - UrpCell{subId} * uPLoc - UrdCell{subId} * lambdaLoc;
                
                pLoc = obj.primalDofsLocal{subId};
                dLoc = obj.dualDofsLocal{subId};
                rLoc = obj.remDofsLocal{subId};
                
                T = obj.Tlocal{subId};
                uLocTransformed = zeros(size(T, 1), 1);
                uLocTransformed(pLoc)         = uPLoc;
                uLocTransformed([rLoc; dLoc]) = uRemLoc;
                
                uLocPhysical = T * uLocTransformed;
                
                uFull(pGlobal)            = uLocPhysical(pLoc);
                uFull([rGlobal; dGlobal]) = uLocPhysical([rLoc; dLoc]);
            end
        end
    end
    
    % =========================================================
    % PUBLIC METHODS: PRECONDITIONER
    % =========================================================
    methods (Access = public)
        
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
                
                Sdd  = obj.localSchurBlocks{subId};
                w    = 1 ./ multiplicity(dualRows);
                rLoc = dualSigns .* r(dualRows);
                
                z(dualRows) = z(dualRows) + dualSigns .* (w .* (Sdd * (w .* rLoc)));
            end
        end
        
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
                
                Sdd = obj.localSchurBlocks{subId};
                w   = 1 ./ multiplicity(dualRows);
                
                W_signed = diag(w .* dualSigns);
                M(dualRows, dualRows) = M(dualRows, dualRows) + W_signed * Sdd * W_signed;
            end
        end
    end
    
    % =========================================================
    % PUBLIC METHODS: VISUALIZATION
    % =========================================================
    methods (Access = public)
        function visualizeFetiNodes(obj)
            allPrimalDofs = unique(vertcat(obj.primalDofsGlobal{:}));
            allDualDofs   = unique(vertcat(obj.dualDofsGlobal{:}));
            allRemDofs    = unique(vertcat(obj.remDofsGlobal{:}));
            
            primalNodes = unique(ceil(allPrimalDofs / obj.dofsPerNode));
            dualNodes   = unique(ceil(allDualDofs   / obj.dofsPerNode));
            remNodes    = unique(ceil(allRemDofs    / obj.dofsPerNode));
            dualNodes   = setdiff(dualNodes, primalNodes);
            remNodes    = setdiff(remNodes,  union(primalNodes, dualNodes));
            
            pCoords = obj.meshCoords(primalNodes, :);
            dCoords = obj.meshCoords(dualNodes,   :);
            rCoords = obj.meshCoords(remNodes,    :);
            
            figure('Name', 'FETI-DP Nodes (Elasticity)', 'Color', 'w');
            hold on; axis equal;
            
            for i = 1:obj.numSubdomains
                locCoords = obj.localMeshes{i}.coord;
                locConnec = obj.localMeshes{i}.connec;                
                patch('Faces', locConnec, 'Vertices', locCoords, ...
                      'FaceColor', 'none', 'EdgeColor', [0.8 0.8 0.8], ...
                      'LineWidth', 1.5, 'HandleVisibility', 'off');
            end
            
            scatter(rCoords(:,1), rCoords(:,2), 20, [0.5 0.5 0.5], 'filled', 'DisplayName', 'Interior');
            scatter(dCoords(:,1), dCoords(:,2), 40, 'b',           'filled', 'DisplayName', 'Interface (Dual)');
            scatter(pCoords(:,1), pCoords(:,2), 80, 'r',           'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1, 'DisplayName', 'Corners (Primal)');
            
            legend('Location', 'bestoutside');
            title('FETI-DP Node Distribution (2D Elasticity)');
            xlabel('X'); ylabel('Y'); 
            grid on; hold off;
        end
    end
    
    % =========================================================
    % PRIVATE METHODS: SUBDOMAIN & MATH OPERATIONS
    % =========================================================
    methods (Access = private)
        
        function extractFetiDofs(obj)
            nSub = obj.numSubdomains;
            tol  = obj.nodeTol;
            
            pGlobal = cell(nSub,1); 
            dGlobal = cell(nSub,1); 
            rGlobal = cell(nSub,1);
            pLocal  = cell(nSub,1); 
            dLocal  = cell(nSub,1); 
            rLocal  = cell(nSub,1);
            eLocal  = cell(nSub,1);
            
            nodeMultiplicity = zeros(size(obj.meshCoords,1), 1);
            for i = 1:nSub
                coords = obj.localMeshes{i}.coord;
                [~, gNodes] = ismembertol(coords, obj.meshCoords, tol, 'ByRows', true);
                nodeMultiplicity(gNodes) = nodeMultiplicity(gNodes) + 1;
            end
            
            isDirNode = false(size(obj.meshCoords,1), 1);
            if ~isempty(obj.dirichletDofs)
                dirNodes = ceil(obj.dirichletDofs / obj.dofsPerNode);
                isDirNode(dirNodes) = true;
            end
            
            for i = 1:nSub
                coords = obj.localMeshes{i}.coord;
                [~, globalNodes] = ismembertol(coords, obj.meshCoords, tol, 'ByRows', true);
                
                minX = min(coords(:,1)); 
                maxX = max(coords(:,1));
                minY = min(coords(:,2)); 
                maxY = max(coords(:,2));
                
                isBL = abs(coords(:,1)-minX) < tol & abs(coords(:,2)-minY) < tol;
                isBR = abs(coords(:,1)-maxX) < tol & abs(coords(:,2)-minY) < tol;
                isTL = abs(coords(:,1)-minX) < tol & abs(coords(:,2)-maxY) < tol;
                isTR = abs(coords(:,1)-maxX) < tol & abs(coords(:,2)-maxY) < tol;
                
                localPrimalNodes = find(isBL | isBR | isTL | isTR);
                
                isShared       = nodeMultiplicity(globalNodes) > 1;
                isDir          = isDirNode(globalNodes);
                localDualNodes = setdiff(find(isShared | isDir), localPrimalNodes);
                localRemNodes  = setdiff(1:size(coords,1), [localPrimalNodes; localDualNodes]);
                
                isOnLeft   = abs(coords(localDualNodes, 1) - minX) < tol;
                isOnRight  = abs(coords(localDualNodes, 1) - maxX) < tol;
                isOnBottom = abs(coords(localDualNodes, 2) - minY) < tol;
                isOnTop    = abs(coords(localDualNodes, 2) - maxY) < tol;
                
                edges = {};
                
                if any(isOnLeft)
                    eNodes = localDualNodes(isOnLeft);
                    [~, sortIdx] = sort(globalNodes(eNodes)); 
                    edges{end+1} = eNodes(sortIdx);
                end
                if any(isOnRight)
                    eNodes = localDualNodes(isOnRight);
                    [~, sortIdx] = sort(globalNodes(eNodes));
                    edges{end+1} = eNodes(sortIdx);
                end
                if any(isOnBottom)
                    eNodes = localDualNodes(isOnBottom);
                    [~, sortIdx] = sort(globalNodes(eNodes));
                    edges{end+1} = eNodes(sortIdx);
                end
                if any(isOnTop)
                    eNodes = localDualNodes(isOnTop);
                    [~, sortIdx] = sort(globalNodes(eNodes));
                    edges{end+1} = eNodes(sortIdx);
                end
                
                eLocal{i} = edges;
                
                masterNodes = [];
                for e = 1:length(edges)
                    masterNodes = [masterNodes; edges{e}(1)]; 
                end
                
                localDualNodes   = setdiff(localDualNodes, masterNodes);
                localPrimalNodes = unique([localPrimalNodes; masterNodes]);
                
                pGlobal{i} = obj.nodesToDofs(globalNodes(localPrimalNodes));
                dGlobal{i} = obj.nodesToDofs(globalNodes(localDualNodes));
                rGlobal{i} = obj.nodesToDofs(globalNodes(localRemNodes));
                
                pLocal{i} = obj.nodesToDofs(localPrimalNodes);
                dLocal{i} = obj.nodesToDofs(localDualNodes);
                rLocal{i} = obj.nodesToDofs(localRemNodes);
            end
            
            obj.primalDofsGlobal = pGlobal; 
            obj.primalDofsLocal  = pLocal;
            obj.dualDofsGlobal   = dGlobal; 
            obj.dualDofsLocal    = dLocal;
            obj.remDofsGlobal    = rGlobal; 
            obj.remDofsLocal     = rLocal;
            obj.edgeNodesLocal   = eLocal; 
        end
        
        function applyEdgeTransformation(obj)
            obj.Tlocal = cell(obj.numSubdomains, 1);
            
            for subId = 1:obj.numSubdomains
                kMat  = obj.localStiffness{subId};
                fVec  = obj.localForces{subId};
                nDofs = length(fVec);
                
                T = speye(nDofs);
                
                edges = obj.edgeNodesLocal{subId};
                
                for e = 1:length(edges)
                    edgeNodes  = edges{e};
                    
                    masterNode = edgeNodes(1);
                    slaveNodes = edgeNodes(2:end);
                    
                    for d = 1:obj.dofsPerNode
                        mDof  = (masterNode - 1) * obj.dofsPerNode + d; 
                        sDofs = (slaveNodes - 1) * obj.dofsPerNode + d; 
                        
                        T(mDof, sDofs) = -1;
                        T(sDofs, mDof) =  1;
                    end
                end
                
                obj.localStiffness{subId} = T' * kMat * T;
                obj.localForces{subId}    = T' * fVec;
                
                obj.Tlocal{subId} = T;
            end
        end
        
        function dofs = nodesToDofs(obj, nodes)
            nodes = nodes(:);
            dofs  = zeros(length(nodes) * obj.dofsPerNode, 1);
            for d = 1:obj.dofsPerNode
                dofs(d:obj.dofsPerNode:end) = (nodes - 1) * obj.dofsPerNode + d;
            end
        end
        
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
        
        function activeIdx = getActivePrimalDofs(obj, allPrimals)
            isDirichlet = ismember(allPrimals, obj.dirichletDofs);
            activeIdx   = find(~isDirichlet);
        end
        
        function Sdd = computeLocalSchur(obj, subId)
            kMat = obj.localStiffness{subId};
            rLoc = obj.remDofsLocal{subId};
            dLoc = obj.dualDofsLocal{subId};
            
            Kii = kMat(rLoc, rLoc);
            Kid = kMat(rLoc, dLoc);
            Kdi = kMat(dLoc, rLoc);
            Kdd = kMat(dLoc, dLoc);
            
            Sdd = Kdd - Kdi * (Kii \ full(Kid));
        end
    end
end