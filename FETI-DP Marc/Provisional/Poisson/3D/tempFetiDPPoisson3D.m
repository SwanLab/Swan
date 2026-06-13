classdef tempFetiDPPoisson3D < handle
    % FETI-DP algebraic solver for 3D Poisson equation on a domain decomposition.
    
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
        
        function obj = tempFetiDPPoisson3D(globalMesh, subMeshes, stiffness, forces, tol, dofsNode)
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
                
                for ii = 1:length(dualRows)
                    for jj = 1:length(dualRows)
                        M(dualRows(ii), dualRows(jj)) = M(dualRows(ii), dualRows(jj)) + ...
                            dualSigns(ii) * w(ii) * Sdd(ii, jj) * w(jj) * dualSigns(jj);
                    end
                end
            end
        end
        
        function visualizeFetiNodes(obj)
            figure('Name', 'FETI-DP Node Classification 3D', 'Color', 'w');
            hold on; axis equal; view(3); grid on;
            
            coords = obj.meshCoords;
            
            for subId = 1:obj.numSubdomains
                localCoords = obj.localMeshes{subId}.coord;
                plot3(localCoords(:, 1), localCoords(:, 2), localCoords(:, 3), 'k.', 'MarkerSize', 1);
            end
            
            for subId = 1:obj.numSubdomains
                pNodes = unique(ceil(obj.primalDofsGlobal{subId} / obj.dofsPerNode));
                dNodes = unique(ceil(obj.dualDofsGlobal{subId} / obj.dofsPerNode));
                rNodes = unique(ceil(obj.remDofsGlobal{subId} / obj.dofsPerNode));
                
                if ~isempty(pNodes)
                    plot3(coords(pNodes, 1), coords(pNodes, 2), coords(pNodes, 3), 'ro', 'MarkerSize', 8, 'LineWidth', 2);
                end
                if ~isempty(dNodes)
                    plot3(coords(dNodes, 1), coords(dNodes, 2), coords(dNodes, 3), 'bs', 'MarkerSize', 5, 'LineWidth', 1.5);
                end
                % Comentado para no saturar el plot en 3D
                % if ~isempty(rNodes)
                %     plot3(coords(rNodes, 1), coords(rNodes, 2), coords(rNodes, 3), 'g^', 'MarkerSize', 3);
                % end
            end
            
            legend('Mesh', 'Primal (corners)', 'Dual (faces/edges)', 'Location', 'best');
            title('FETI-DP DOF Classification for 3D Poisson');
            xlabel('X'); ylabel('Y'); zlabel('Z');
            hold off;
        end

        % function extractFetiDofs(obj)
        %     nSub = prod(obj.numSubdomains);
        %     tol  = obj.nodeTol;
        % 
        %     pGlobal = cell(nSub,1);
        %     dGlobal = cell(nSub,1);
        %     rGlobal = cell(nSub,1);
        % 
        %     pLocal = cell(nSub,1);
        %     dLocal = cell(nSub,1);
        %     rLocal = cell(nSub,1);
        % 
        %     for i = 1:nSub
        %         coords  = obj.localMeshes{i}.coord;
        %         nNodes  = size(coords,1);
        % 
        %         % Local -> global map
        %         [~,globalNodes] = ismembertol(coords, obj.meshCoords, tol, 'ByRows', true);
        % 
        %         % Boundary meshes (en 3D devolverá la malla de superficie externa del volumen)
        %         bMeshes = obj.localMeshes{i}.createBoundaryMesh();
        %         gcCell = cell(length(bMeshes), 1);
        %         for b = 1:length(bMeshes)
        %             gcCell{b} = bMeshes{b}.globalConnec(:);
        %         end
        %         boundaryNodes = unique(vertcat(gcCell{:}));
        % 
        %         % --- EXTRACCIÓN ROBUSTA DE ESQUINAS (PRIMALES) PARA 3D ---
        %         minX = min(coords(:, 1)); maxX = max(coords(:, 1));
        %         minY = min(coords(:, 2)); maxY = max(coords(:, 2));
        %         minZ = min(coords(:, 3)); maxZ = max(coords(:, 3));
        % 
        %         % Una esquina en 3D es el punto de intersección de 3 planos extremos (X, Y y Z)
        %         isCorner = (abs(coords(:, 1) - minX) < tol | abs(coords(:, 1) - maxX) < tol) & ...
        %                    (abs(coords(:, 2) - minY) < tol | abs(coords(:, 2) - maxY) < tol) & ...
        %                    (abs(coords(:, 3) - minZ) < tol | abs(coords(:, 3) - maxZ) < tol);
        % 
        %         primalNodes = find(isCorner);
        % 
        %         % Dual = boundary minus primal
        %         dualNodes = setdiff(boundaryNodes, primalNodes);
        % 
        %         % Remaining = internal
        %         remNodes = setdiff((1:nNodes)', [primalNodes; dualNodes]);
        % 
        %         % Store local
        %         pLocal{i} = primalNodes;
        %         dLocal{i} = dualNodes;
        %         rLocal{i} = remNodes;
        % 
        %         % Store global
        %         pGlobal{i} = globalNodes(primalNodes);
        %         dGlobal{i} = globalNodes(dualNodes);
        %         rGlobal{i} = globalNodes(remNodes);
        %     end
        % 
        %     obj.primalDofsGlobal = pGlobal;
        %     obj.dualDofsGlobal   = dGlobal;
        %     obj.remDofsGlobal    = rGlobal;
        % 
        %     obj.primalDofsLocal  = pLocal;
        %     obj.dualDofsLocal    = dLocal;
        %     obj.remDofsLocal     = rLocal;
        % end
        
        function extractFetiDofs(obj)
            nSub = prod(obj.numSubdomains);
            tol  = obj.nodeTol;
            
            pGlobal = cell(nSub,1);
            dGlobal = cell(nSub,1);
            rGlobal = cell(nSub,1);
            
            pLocal = cell(nSub,1);
            dLocal = cell(nSub,1);
            rLocal = cell(nSub,1);
            
            for i = 1:nSub
                coords  = obj.localMeshes{i}.coord;
                nNodes  = size(coords,1);
                
                % Mapeo Local -> Global
                [~,globalNodes] = ismembertol(coords, obj.meshCoords, tol, 'ByRows', true);
                
                % Calculamos los límites extremos del subdominio actual
                minX = min(coords(:, 1)); maxX = max(coords(:, 1));
                minY = min(coords(:, 2)); maxY = max(coords(:, 2));
                minZ = min(coords(:, 3)); maxZ = max(coords(:, 3));
                
                % Detectamos si un nodo pertenece a un plano exterior
                onX = (abs(coords(:, 1) - minX) < tol | abs(coords(:, 1) - maxX) < tol);
                onY = (abs(coords(:, 2) - minY) < tol | abs(coords(:, 2) - maxY) < tol);
                onZ = (abs(coords(:, 3) - minZ) < tol | abs(coords(:, 3) - maxZ) < tol);
                
                % Sumamos a cuántos planos extremos pertenece cada nodo
                % 1 = Cara, 2 = Arista, 3 = Esquina
                numFaces = onX + onY + onZ;
                
                % Frontera total (cualquier nodo que toque al menos 1 plano)
                boundaryNodes = find(numFaces >= 1);
                
                % --- SOLUCIÓN AL BUG 3D (Wirebasket) ---
                % Primales = Esquinas (3) + Aristas (2)
                primalNodes = find(numFaces >= 3);
                
                % Duales = Solo nodos interiores a las caras (exactamente 1)
                % Esto garantiza que siempre sean compartidos por solo 2 subdominios.
                dualNodes = find(numFaces >= 1 & numFaces < 3);
                
                % Restantes = Nodos interiores al volumen
                remNodes = setdiff((1:nNodes)', boundaryNodes);
                
                % Almacenamos local
                pLocal{i} = primalNodes;
                dLocal{i} = dualNodes;
                rLocal{i} = remNodes;
                
                % Almacenamos global
                pGlobal{i} = globalNodes(primalNodes);
                dGlobal{i} = globalNodes(dualNodes);
                rGlobal{i} = globalNodes(remNodes);
            end
            
            obj.primalDofsGlobal = pGlobal;
            obj.dualDofsGlobal   = dGlobal;
            obj.remDofsGlobal    = rGlobal;
            
            obj.primalDofsLocal  = pLocal;
            obj.dualDofsLocal    = dLocal;
            obj.remDofsLocal     = rLocal;
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