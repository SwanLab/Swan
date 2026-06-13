classdef FetidpV2 < handle
    % FETIDP Classe encarregada de la lògica algebraica del mètode FETI-DP.
    
    properties (Access = private)
        localStiffness  
        localForces    
        localMeshes   
        meshCoords      
        dofsPerNode     
        nodeTol         
        numSubdomains   
        dirichletDofs  
        
        primalDofsGlobal      
        dualDofsGlobal        
        remDofsGlobal         
        
        primalDofsLocal      
        dualDofsLocal        
        remDofsLocal         
    end
    
    properties (Access = public)
        bdLocal         
    end
    
    methods (Access = public)
        function obj = FetidpV2(globalMesh, subMeshes, stiffness, forces, tol, dofsNode, dirDofs)
            obj.meshCoords     = globalMesh.coord;
            obj.dofsPerNode    = dofsNode;
            obj.localMeshes    = subMeshes;
            obj.localStiffness = stiffness;
            obj.localForces    = forces;
            obj.nodeTol        = tol;
            obj.numSubdomains  = length(stiffness);
            obj.dirichletDofs  = dirDofs;
            obj.bdLocal        = cell(obj.numSubdomains, 1);
            
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
                
                Urd = Krr \ full(Tdr');     
                Urp = Krr \ full(Krp);             
                
                ap = obj.createAp(subId, allPrimals);
                [bd, visitedDuals] = obj.createBd(subId, allDuals, visitedDuals);
                obj.bdLocal{subId} = bd;
                                
                SppLoc = Kpp - Kpr * Urp;
                SPP = SPP + ap * SppLoc * ap';
                
                if size(bd, 2) > 0
                    fDual = fDual + bd * (Urd' * Tdr') * bd';
                    dBar  = dBar + bd * (Urd' * fR);
                    BrKrrInvKrp = BrKrrInvKrp + bd * (Tdr * Urp) * ap';
                end
                
                rhsPrimal = rhsPrimal + ap * (fP - Urp' * fR);
            end
            
            activeDofs = obj.getActivePrimalDofs(allPrimals);   
            
            SPPActive = SPP(activeDofs, activeDofs);                    
            BrKrrInvKrpActive  = BrKrrInvKrp(:, activeDofs);               
            rhsActive = rhsPrimal(activeDofs);                     
            
            fMat = fDual + BrKrrInvKrpActive * (SPPActive \ full(BrKrrInvKrpActive'));    
            dBar = dBar  - BrKrrInvKrpActive * (SPPActive \ rhsActive);         
        end
        
        % function zVec = applyDirichletPrecond(obj, rVec)
        %     numL = length(rVec);
        %     zVec = zeros(numL, 1);
        %     scaling = 0.5; 
        % 
        %     for subId = 1:obj.numSubdomains
        %         bd = obj.bdLocal{subId};
        %         if size(bd, 2) == 0
        %             continue; 
        %         end
        % 
        %         kMat = obj.localStiffness{subId};
        %         dLoc = obj.dualDofsLocal{subId};
        %         rLoc = obj.remDofsLocal{subId}; 
        % 
        %         kIi = kMat(rLoc, rLoc);
        %         kId = kMat(rLoc, dLoc);
        %         kDi = kMat(dLoc, rLoc);
        %         kDd = kMat(dLoc, dLoc);
        % 
        %         sLocal = kDd - kDi * (kIi \ kId);
        %         rLocal = scaling * (bd' * rVec);
        %         zLocal = sLocal * rLocal;
        % 
        %         zVec = zVec + bd * (scaling * zLocal);
        %     end
        % end
        
        function uFull = reconstructGlobalSolution(obj, lambdaSol, numNodes)
            allPrimals = unique(vertcat(obj.primalDofsGlobal{:}));
            numP = length(allPrimals);
            SPP = sparse(numP, numP);
            rhsPrimal = zeros(numP, 1);
            
            UrpCell = cell(obj.numSubdomains, 1);
            UrdCell = cell(obj.numSubdomains, 1);
            KrrInvFrCell = cell(obj.numSubdomains, 1);
            
            for subId = 1:obj.numSubdomains
                [Krr, Krp, Kpr, Kpp, fR, fP, Tdr] = obj.splitLocalMatrices(subId);
                ap = obj.createAp(subId, allPrimals);
                bd = obj.bdLocal{subId};
                
                Urp = Krr \ full(Krp); 
                Urd = Krr \ full(Tdr');
                KrrInvFr = Krr \ fR;
                
                UrpCell{subId} = Urp;
                UrdCell{subId} = Urd;
                KrrInvFrCell{subId} = KrrInvFr;
                
                sppLoc = Kpp - Kpr * Urp;
                SPP = SPP + ap * sppLoc * ap';
                
                term = Tdr' * (bd' * lambdaSol) - fR;
                rhsPrimal = rhsPrimal + ap * (fP + Urp' * term);
            end
            
            activeDofs = obj.getActivePrimalDofs(allPrimals);
            sppActive = SPP(activeDofs, activeDofs);
            rhsActive = rhsPrimal(activeDofs);
            
            uP = zeros(numP, 1);
            uP(activeDofs) = sppActive \ rhsActive;
            
            uFull = zeros(numNodes * obj.dofsPerNode, 1);
            
            for subId = 1:obj.numSubdomains
                pGlobal = obj.primalDofsGlobal{subId};
                dGlobal = obj.dualDofsGlobal{subId};
                rGlobal = obj.remDofsGlobal{subId};
                
                ap = obj.createAp(subId, allPrimals);
                bd = obj.bdLocal{subId};
                
                Urp = UrpCell{subId};
                Urd = UrdCell{subId};
                KrrInvFr = KrrInvFrCell{subId};
                
                uPLoc = ap' * uP;
                uRemLoc = KrrInvFr - Urp * uPLoc - Urd * (bd' * lambdaSol);
                
                uFull(pGlobal) = uPLoc;
                
                allRemGlobalSorted = [rGlobal; dGlobal];
                uFull(allRemGlobalSorted) = uRemLoc;
            end
        end
        
        function visualizeFetiNodes(obj)

            allPrimalDofs = unique(vertcat(obj.primalDofsGlobal{:}));
            allDualDofs   = unique(vertcat(obj.dualDofsGlobal{:}));
            allRemDofs    = unique(vertcat(obj.remDofsGlobal{:}));
            
            primalNodes = unique(ceil(allPrimalDofs / obj.dofsPerNode));
            dualNodes   = unique(ceil(allDualDofs / obj.dofsPerNode));
            remNodes    = unique(ceil(allRemDofs / obj.dofsPerNode));
            
            dualNodes = setdiff(dualNodes, primalNodes);
            remNodes  = setdiff(remNodes, union(primalNodes, dualNodes));
            
            pCoords = obj.meshCoords(primalNodes, :);
            dCoords = obj.meshCoords(dualNodes, :);
            rCoords = obj.meshCoords(remNodes, :);
            
            figure('Name', 'Nodes FETI-DP', 'Color', 'w');
            hold on; axis equal;
            
            scatter(rCoords(:,1), rCoords(:,2), 20, [0.5 0.5 0.5], 'filled', 'DisplayName', 'Interns (Remaining)');
            scatter(dCoords(:,1), dCoords(:,2), 40, 'b', 'filled', 'DisplayName', 'Interfície (Dual)');
            scatter(pCoords(:,1), pCoords(:,2), 80, 'r', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1, 'DisplayName', 'Cantonades (Primal)');
            
            legend('Location', 'bestoutside');
            title('Distribució de Nodes FETI-DP');
            xlabel('X'); ylabel('Y');
            grid on;
            hold off;
        end

        function z = applyDirichletPrecond(obj, r)
            allPrimals = unique(vertcat(obj.primalDofsGlobal{:}));
            allDuals   = setdiff(unique(vertcat(obj.dualDofsGlobal{:})), allPrimals);
            numDuals   = length(allDuals);
            
            % Calcula multiplicitat de cada dual global
            multiplicity = zeros(numDuals, 1);
            for subId = 1:obj.numSubdomains
                dualRows = obj.dualIdxLocal{subId};
                multiplicity(dualRows) = multiplicity(dualRows) + 1;
            end

            z = zeros(numDuals, 1);
        
            for subId = 1:obj.numSubdomains
                Sdd = obj.computeLocalSchur(subId);   
        
                dualRows  = obj.dualIdxLocal{subId};
                dualSigns = obj.dualSignsLocal{subId};
        
                if isempty(dualRows), continue; end

                % Pesos locals
                w = 1./multiplicity(dualRows);
        
                % Extreu la part local de r (amb signe de Bd^T)
                rLoc = dualSigns .* r(dualRows);
        
                % Aplica el Schur local i escampa (amb signe de Bd)
                z(dualRows) = z(dualRows) + dualSigns .* (w .* (Sdd *(w .* rLoc)));
            end
        end
    end
    
    methods (Access = private)
        function extractFetiDofs(obj)
            nSub = prod(obj.numSubdomains);
            tol = obj.nodeTol;
            
            pGlobal = cell(nSub, 1); 
            dGlobal = cell(nSub, 1); 
            rGlobal = cell(nSub, 1);

            pLocal = cell(nSub, 1);  
            dLocal = cell(nSub, 1);  
            rLocal = cell(nSub, 1);
            
            nodeMultiplicity = zeros(size(obj.meshCoords, 1), 1);
            for i = 1:nSub
                coords = obj.localMeshes{i}.coord;
                [~, globalNodes] = ismembertol(coords, obj.meshCoords, tol, 'ByRows', true);
                nodeMultiplicity(globalNodes) = nodeMultiplicity(globalNodes) + 1;
            end
            
            isDirNode = false(size(obj.meshCoords, 1), 1);
            if ~isempty(obj.dirichletDofs)
                dirNodes = ceil(obj.dirichletDofs / obj.dofsPerNode);
                isDirNode(dirNodes) = true;
            end
            
            for i = 1:nSub
                coords = obj.localMeshes{i}.coord;
                [~, globalNodes] = ismembertol(coords, obj.meshCoords, tol, 'ByRows', true);
                
                minX = min(coords(:,1)); maxX = max(coords(:,1));
                minY = min(coords(:,2)); maxY = max(coords(:,2));
                
                isBL = abs(coords(:,1) - minX) < tol & abs(coords(:,2) - minY) < tol;
                isBR = abs(coords(:,1) - maxX) < tol & abs(coords(:,2) - minY) < tol;
                isTL = abs(coords(:,1) - minX) < tol & abs(coords(:,2) - maxY) < tol;
                isTR = abs(coords(:,1) - maxX) < tol & abs(coords(:,2) - maxY) < tol;
                
                localPrimalNodes = find(isBL | isBR | isTL | isTR);
                
                isShared = nodeMultiplicity(globalNodes) > 1;
                isDir    = isDirNode(globalNodes);
                isDual   = isShared | isDir;
                
                localDualNodes = setdiff(find(isDual), localPrimalNodes);
                localRemNodes  = setdiff(1:size(coords,1), [localPrimalNodes; localDualNodes]);
                
                % GUARDEM ELS DOFS GLOBALS I LOCALS EXACTES (Sense usar ismember desprès)
                pGlobal{i} = obj.nodesToDofs(globalNodes(localPrimalNodes));
                dGlobal{i} = obj.nodesToDofs(globalNodes(localDualNodes));
                rGlobal{i} = obj.nodesToDofs(globalNodes(localRemNodes));
                
                pLocal{i} = obj.nodesToDofs(localPrimalNodes);
                dLocal{i} = obj.nodesToDofs(localDualNodes);
                rLocal{i} = obj.nodesToDofs(localRemNodes);
            end
            
            obj.primalDofsGlobal = pGlobal; obj.primalDofsLocal = pLocal;
            obj.dualDofsGlobal   = dGlobal; obj.dualDofsLocal   = dLocal;
            obj.remDofsGlobal    = rGlobal; obj.remDofsLocal    = rLocal;          
        end
        
        function dofs = nodesToDofs(obj, nodes)            
            nodes = nodes(:);
            dofs = zeros(length(nodes) * obj.dofsPerNode, 1);
            for d = 1:obj.dofsPerNode
                dofs(d:obj.dofsPerNode:end) = (nodes-1) * obj.dofsPerNode + d;
            end
        end
        
        function [Krr, Krp, Kpr, Kpp, fR, fP, Tdr] = splitLocalMatrices(obj, subId)
            kMat = obj.localStiffness{subId};
            fVec = obj.localForces{subId};
            
            pLoc = obj.primalDofsLocal{subId};
            dLoc = obj.dualDofsLocal{subId};
            rLoc = obj.remDofsLocal{subId};
            
            % Ajuntem els DOFs sense desordenar-los (evitant el sort)
            remLoc = [rLoc; dLoc]; 
            
            Krr = kMat(remLoc, remLoc);
            Krp = kMat(remLoc, pLoc);
            Kpr = kMat(pLoc, remLoc); 
            Kpp = kMat(pLoc, pLoc);
            
            fR = fVec(remLoc);
            fP = fVec(pLoc);
                        
            numD = length(dLoc); 
            numR = length(remLoc);      
            
            rows = (1:numD)';
            cols = (length(rLoc) + 1 : numR)';
            vals = ones(numD, 1);
            
            Tdr = sparse(rows, cols, vals, numD, numR);
        end
        
        function ap = createAp(obj, subId, allPrimals)
            pGlobal = obj.primalDofsGlobal{subId};
            numPrimalsGlobal = length(allPrimals);
            numPrimalsLocal  = length(pGlobal);
            
            [~, rows] = ismember(pGlobal, allPrimals);
            cols = (1:numPrimalsLocal)';
            vals = ones(numPrimalsLocal, 1);
            
            ap = sparse(rows, cols, vals, numPrimalsGlobal, numPrimalsLocal);
        end
        
        function [bd, visitedDuals] = createBd(obj, subId, allDuals, visitedDuals)
            dGlobal = obj.dualDofsGlobal{subId};
            numDualsGlobal = length(allDuals);
            numDualsLocal  = length(dGlobal);
        
            if numDualsLocal == 0
                bd = sparse(numDualsGlobal, 0);
                return;
            end
        
            [~, rows] = ismember(dGlobal, allDuals);
            cols = (1:numDualsLocal)';
        
            isFirstVisit = (visitedDuals(dGlobal) == 0);
            vals = ones(numDualsLocal, 1);
            vals(~isFirstVisit) = -1;
        
            visitedDuals(dGlobal(isFirstVisit)) = 1;
        
            bd = sparse(rows, cols, vals, numDualsGlobal, numDualsLocal);
        end
        
        function activeIdx = getActivePrimalDofs(obj, allPrimals)
            isDirichlet = ismember(allPrimals, obj.dirichletDofs);
            activeIdx = find(~isDirichlet);  
        end

        function Sdd = computeLocalSchur(obj, subId)
            kMat = obj.localStiffness{subId};
            rLoc = obj.remDofsLocal{subId};   % interiors
            dLoc = obj.dualDofsLocal{subId};  % interfície dual
        
            Kii = kMat(rLoc, rLoc);
            Kid = kMat(rLoc, dLoc);
            Kdi = kMat(dLoc, rLoc);
            Kdd = kMat(dLoc, dLoc);
        
            % S_dd = Kdd - Kdi * Kii^{-1} * Kid
            Sdd = Kdd - Kdi * (Kii \ full(Kid));
        end
    end
end