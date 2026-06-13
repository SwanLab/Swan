classdef tempFetiDPPoisson < handle
    % FETI-DP algebraic solver for 2D Poisson equation on a domain decomposition.
    
    properties (Access = private)
        localStiffness
        localForces
        localMeshes
        meshCoords
        dofsPerNode
        nodeTol
        numSubdomains
        
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
        
        function obj = tempFetiDPPoisson(globalMesh, subMeshes, stiffness, forces, tol, dofsNode)
            obj.meshCoords          = globalMesh.coord;
            obj.dofsPerNode         = dofsNode;
            obj.localMeshes         = subMeshes;
            obj.localStiffness      = stiffness;
            obj.localForces         = forces;
            obj.nodeTol             = tol;
            obj.numSubdomains       = length(stiffness);
            obj.primalIdxLocal      = cell(obj.numSubdomains, 1);
            
            obj.extractFetiDofs();
        end
        
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
        
        function M = buildPrecondMatrix(obj)
            allPrimals = unique(vertcat(obj.primalDofsGlobal{:}));
            allDuals   = setdiff(unique(vertcat(obj.dualDofsGlobal{:})), allPrimals);
            numDuals   = length(allDuals);
            
            multiplicity = zeros(numDuals, 1);
            for subId = 1:obj.numSubdomains
                dualRows = obj.dualIdxLocal{subId};
                multiplicity(dualRows) = multiplicity(dualRows) + 1;
            end
            
            M = sparse(numDuals, numDuals);
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
        
        function visualizeFetiNodes(obj)
            figure('Name', 'FETI-DP Node Classification (Poisson)', 'Color', 'w');
            hold on; axis equal;
            
            coords = obj.meshCoords;
            
            for subId = 1:obj.numSubdomains
                localCoords = obj.localMeshes{subId}.coord;
                plot(localCoords(:, 1), localCoords(:, 2), 'k.', 'MarkerSize', 1);
            end
            
            for subId = 1:obj.numSubdomains
                pNodes = unique(ceil(obj.primalDofsGlobal{subId} / obj.dofsPerNode));
                dNodes = unique(ceil(obj.dualDofsGlobal{subId} / obj.dofsPerNode));
                rNodes = unique(ceil(obj.remDofsGlobal{subId} / obj.dofsPerNode));
                
                if ~isempty(pNodes)
                    plot(coords(pNodes, 1), coords(pNodes, 2), 'ro', 'MarkerSize', 8, 'LineWidth', 2);
                end
                if ~isempty(dNodes)
                    plot(coords(dNodes, 1), coords(dNodes, 2), 'bs', 'MarkerSize', 6, 'LineWidth', 1.5);
                end
                if ~isempty(rNodes)
                    plot(coords(rNodes, 1), coords(rNodes, 2), 'g^', 'MarkerSize', 4);
                end
            end
            
            legend('Mesh', 'Primal (corners)', 'Dual (exterior edges)', 'Remaining (interior)', 'Location', 'best');
            title('FETI-DP DOF Classification for Poisson');
            xlabel('X'); ylabel('Y');
            hold off;
        end

        function extractFetiDofs(obj)
            nSub = obj.numSubdomains; 
            tol  = obj.nodeTol;
        
            pGlobal = cell(nSub,1);
            dGlobal = cell(nSub,1);
            rGlobal = cell(nSub,1);
        
            pLocal = cell(nSub,1);
            dLocal = cell(nSub,1);
            rLocal = cell(nSub,1);
        
            for i = 1:nSub
                coords = obj.localMeshes{i}.coord;
                nNodes = size(coords,1);
        
                [~, globalNodes] = ismembertol(coords, obj.meshCoords, tol, 'ByRows', true);
        
                bMeshes = obj.localMeshes{i}.createBoundaryMesh();
                boundaryNodes = [];
                primalNodes   = [];
        
                for b = 1:length(bMeshes)
                    gc = bMeshes{b}.globalConnec;
                    boundaryNodes = [boundaryNodes; gc(:)];
        
                    uN  = unique(gc(:));
                    deg = accumarray(gc(:), 1);
                    corners = uN(deg(uN) == 1);
        
                    primalNodes = [primalNodes; corners];
                end
        
                boundaryNodes = unique(boundaryNodes);
                primalNodes   = unique(primalNodes);
        
                % Dual 
                dualNodes = setdiff(boundaryNodes, primalNodes);
        
                % Remaining 
                remNodes = setdiff((1:nNodes)', [primalNodes; dualNodes]);
        
                obj.primalDofsLocal{i} = obj.nodesToDofs(primalNodes);
                obj.dualDofsLocal{i}   = obj.nodesToDofs(dualNodes);
                obj.remDofsLocal{i}    = obj.nodesToDofs(remNodes);
        
                obj.primalDofsGlobal{i} = obj.nodesToDofs(globalNodes(primalNodes));
                obj.dualDofsGlobal{i}   = obj.nodesToDofs(globalNodes(dualNodes));
                obj.remDofsGlobal{i}    = obj.nodesToDofs(globalNodes(remNodes));
            end
        end

           

        function dofs = nodesToDofs(obj, nodes)
            nodes = nodes(:);
            base  = (nodes - 1) * obj.dofsPerNode;
            dofs  = reshape(base + (1:obj.dofsPerNode), [], 1);
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
        
        function pRows = getPrimalRows(obj, subId, allPrimals)
            pRows = obj.primalIdxLocal{subId};
            
            if isempty(pRows)
                [~, pRows] = ismember(obj.primalDofsGlobal{subId}, allPrimals);
            end
        end
        
        function activeIdx = getActivePrimalDofs(~, allPrimals)
            activeIdx = (1:length(allPrimals));
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
    end
end
