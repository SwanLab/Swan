classdef FetiDPPoissonV1 < handle
    % FETIDPPOISSON Classe per a la lògica algebraica de FETI-DP per al problema de Poisson.
    % Replica exactament la lògica de Python: tot el contorn exterior és Dirichlet homogeni (u=0).
    
    properties (Access = private)
        localStiffness  
        localForces     
        localMeshes     
        meshCoords      
        nodeTol         
        numSubdomains   
        
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
        function obj = FetiDPPoissonV1(globalMesh, subMeshes, stiffness, forces, tol)
            obj.meshCoords     = globalMesh.coord;
            obj.localMeshes    = subMeshes;
            obj.localStiffness = stiffness;
            obj.localForces    = forces;
            obj.nodeTol        = tol;
            obj.numSubdomains  = length(stiffness);
            obj.bdLocal        = cell(obj.numSubdomains, 1);
            
            % Classifiquem els graus de llibertat com a Python
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
                        
            warning('off', 'MATLAB:singularMatrix');
            warning('off', 'MATLAB:nearlySingularMatrix');
            
            for subId = 1:obj.numSubdomains
                [Krr, Krp, Kpr, Kpp, fR, fP, Tdr] = obj.splitLocalMatrices(subId);
                
                Urd = Krr \ full(Tdr');     
                Urp = Krr \ full(Krp);             
                
                Ap = obj.createAp(subId, allPrimals);
                [Bd, visitedDuals] = obj.createBd(subId, allDuals, visitedDuals);
                obj.bdLocal{subId} = Bd;
                                
                sppLoc = Kpp - Kpr * Urp;
                SPP = SPP + Ap * sppLoc * Ap';
                
                if size(Bd, 2) > 0
                    fDual = fDual + Bd * (Urd' * Tdr') * Bd';
                    dBar  = dBar + Bd * (Urd' * fR);
                    BrKrrInvKrp = BrKrrInvKrp + Bd * (Tdr * Urp) * Ap';
                end
                
                rhsPrimal = rhsPrimal + Ap * (fP - Urp' * fR);
            end
            
            warning('on', 'MATLAB:singularMatrix');
            warning('on', 'MATLAB:nearlySingularMatrix');
                        
            activeDofs = obj.getActivePrimalDofs(allPrimals);   
            
            SPPActive = SPP(activeDofs, activeDofs);                    
            BrActive  = BrKrrInvKrp(:, activeDofs);               
            rhsActive = rhsPrimal(activeDofs);                     
            
            fMat = fDual + BrActive * (SPPActive \ full(BrActive'));    
            dBar = dBar  - BrActive * (SPPActive \ rhsActive);         
        end
        
        function zVec = applyDirichletPrecond(obj, rVec)
            numL = length(rVec);
            zVec = zeros(numL, 1);
            
            % Calcular la multiplicitat de cada DOF dual per l'escalat
            allDuals = unique(vertcat(obj.dualDofsGlobal{:}));
            multiplicity = zeros(numL, 1);
            for subId = 1:obj.numSubdomains
                dGlobal = obj.dualDofsGlobal{subId};
                [~, posInGlobal] = ismember(dGlobal, allDuals);
                validIdx = posInGlobal > 0;
                multiplicity(posInGlobal(validIdx)) = multiplicity(posInGlobal(validIdx)) + 1;
            end
            
            W = spdiags(1 ./ multiplicity, 0, numL, numL);
            scaledR = W * rVec;
            
            for subId = 1:obj.numSubdomains
                bd = obj.bdLocal{subId};
                if size(bd, 2) == 0
                    continue; 
                end
                
                kMat = obj.localStiffness{subId};
                
                dLoc = obj.dualDofsLocal{subId};
                rLoc = obj.remDofsLocal{subId}; 
                
                kIi = kMat(rLoc, rLoc);
                kId = kMat(rLoc, dLoc);
                kDi = kMat(dLoc, rLoc);
                kDd = kMat(dLoc, dLoc);
                
                sLocal = kDd - kDi * (kIi \ kId);
                
                rLocal = bd' * scaledR;
                zLocal = sLocal * rLocal;
                
                zVec = zVec + W * (bd * zLocal);
            end
        end
        
        function uFull = reconstructGlobalSolution(obj, lambdaSol, numNodes)
            allPrimals = unique(vertcat(obj.primalDofsGlobal{:}));
            numP = length(allPrimals);
            SPP = sparse(numP, numP);
            rhsPrimal = zeros(numP, 1);
            
            UrpCell = cell(obj.numSubdomains, 1);
            UrdCell = cell(obj.numSubdomains, 1);
            KrrInvFrCell = cell(obj.numSubdomains, 1);
                        
            warning('off', 'MATLAB:singularMatrix');
            warning('off', 'MATLAB:nearlySingularMatrix');
            
            for subId = 1:obj.numSubdomains
                [Krr, Krp, Kpr, Kpp, fR, fP, Tdr] = obj.splitLocalMatrices(subId);
                Ap = obj.createAp(subId, allPrimals);
                Bd = obj.bdLocal{subId};
                
                Urp = Krr \ full(Krp); 
                Urd = Krr \ full(Tdr');
                KrrInvFr = Krr \ fR;
                
                UrpCell{subId} = Urp;
                UrdCell{subId} = Urd;
                KrrInvFrCell{subId} = KrrInvFr;
                
                sppLoc = Kpp - Kpr * Urp;
                SPP = SPP + Ap * sppLoc * Ap';
                
                term = Tdr' * (Bd' * lambdaSol) - fR;
                rhsPrimal = rhsPrimal + Ap * (fP + Urp' * term);
            end
            
            warning('on', 'MATLAB:singularMatrix');
            warning('on', 'MATLAB:nearlySingularMatrix');
                        
            activeDofs = obj.getActivePrimalDofs(allPrimals);
            sppActive = SPP(activeDofs, activeDofs);
            rhsActive = rhsPrimal(activeDofs);
            
            uP = zeros(numP, 1);
            uP(activeDofs) = sppActive \ rhsActive;
            
            uFull = zeros(numNodes * 1, 1); 
            
            for subId = 1:obj.numSubdomains
                pGlobal = obj.primalDofsGlobal{subId};
                dGlobal = obj.dualDofsGlobal{subId};
                rGlobal = obj.remDofsGlobal{subId};
                
                Ap = obj.createAp(subId, allPrimals);
                Bd = obj.bdLocal{subId};
                
                Urp = UrpCell{subId};
                Urd = UrdCell{subId};
                KrrInvFr = KrrInvFrCell{subId};
                
                uPLoc = Ap' * uP;
                uRemLoc = KrrInvFr - Urp * uPLoc - Urd * (Bd' * lambdaSol);
                
                uFull(pGlobal) = uPLoc;
                
                % SENSE SORT: Mantenim la mateixa estructura [r; d] que a splitLocalMatrices
                allRemGlobal = [rGlobal; dGlobal];
                uFull(allRemGlobal) = uRemLoc;
            end
        end
        
        function visualizeFetiNodes(obj)
            allPrimalDofs = unique(vertcat(obj.primalDofsGlobal{:}));
            allDualDofs   = unique(vertcat(obj.dualDofsGlobal{:}));
            allRemDofs    = unique(vertcat(obj.remDofsGlobal{:}));
            
            primalNodes = allPrimalDofs; 
            dualNodes   = allDualDofs;
            remNodes    = allRemDofs;
            
            dualNodes = setdiff(dualNodes, primalNodes);
            remNodes  = setdiff(remNodes, union(primalNodes, dualNodes));
            
            pCoords = obj.meshCoords(primalNodes, :);
            dCoords = obj.meshCoords(dualNodes, :);
            rCoords = obj.meshCoords(remNodes, :);
            
            figure('Name', 'Nodes FETI-DP (Poisson)', 'Color', 'w');
            hold on; axis equal;
            
            scatter(rCoords(:,1), rCoords(:,2), 20, [0.5 0.5 0.5], 'filled', 'DisplayName', 'Interns (Remaining)');
            scatter(dCoords(:,1), dCoords(:,2), 40, 'b', 'filled', 'DisplayName', 'Interfície + Contorn (Dual)');
            scatter(pCoords(:,1), pCoords(:,2), 80, 'r', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1, 'DisplayName', 'Cantonades (Primal)');
            
            legend('Location', 'bestoutside');
            title('Distribució de Nodes FETI-DP (Poisson)');
            xlabel('X'); ylabel('Y');
            grid on;
            hold off;
        end
    end
    
    methods (Access = private)
        function extractFetiDofs(obj)
            nSub = prod(obj.numSubdomains);
            tol = obj.nodeTol;
            
            pGlobal = cell(nSub, 1); dGlobal = cell(nSub, 1); rGlobal = cell(nSub, 1);
            pLocal  = cell(nSub, 1); dLocal  = cell(nSub, 1); rLocal  = cell(nSub, 1);
            
            for i = 1:nSub
                coords = obj.localMeshes{i}.coord;
                
                minX = min(coords(:,1)); maxX = max(coords(:,1));
                minY = min(coords(:,2)); maxY = max(coords(:,2));
                
                isBL = abs(coords(:,1) - minX) < tol & abs(coords(:,2) - minY) < tol;
                isBR = abs(coords(:,1) - maxX) < tol & abs(coords(:,2) - minY) < tol;
                isTL = abs(coords(:,1) - minX) < tol & abs(coords(:,2) - maxY) < tol;
                isTR = abs(coords(:,1) - maxX) < tol & abs(coords(:,2) - maxY) < tol;
                
                localPrimalNodes = find(isBL | isBR | isTL | isTR);
                
                % Tota frontera és Dual (Dirichlet global implícit)
                isOnLocalBoundary = abs(coords(:,1) - minX) < tol | abs(coords(:,1) - maxX) < tol | ...
                                    abs(coords(:,2) - minY) < tol | abs(coords(:,2) - maxY) < tol;
                
                isDual = isOnLocalBoundary; 
                localDualNodes = setdiff(find(isDual), localPrimalNodes);
                
                isRem = ~isDual;
                localRemNodes = setdiff(find(isRem), localPrimalNodes);
                
                [~, globalNodes] = ismembertol(coords, obj.meshCoords, tol, 'ByRows', true);                
                
                pGlobal{i} = globalNodes(localPrimalNodes);
                dGlobal{i} = globalNodes(localDualNodes);
                rGlobal{i} = globalNodes(localRemNodes);
                
                pLocal{i} = localPrimalNodes;
                dLocal{i} = localDualNodes;
                rLocal{i} = localRemNodes;
            end
            
            obj.primalDofsGlobal = pGlobal; obj.primalDofsLocal = pLocal;
            obj.dualDofsGlobal   = dGlobal; obj.dualDofsLocal   = dLocal;
            obj.remDofsGlobal    = rGlobal; obj.remDofsLocal    = rLocal;          
        end
        
        function [Krr, Krp, Kpr, Kpp, fR, fP, Tdr] = splitLocalMatrices(obj, subId)
            kMat = obj.localStiffness{subId};
            fVec = obj.localForces{subId};
            
            pLoc = obj.primalDofsLocal{subId};
            dLoc = obj.dualDofsLocal{subId};
            rLoc = obj.remDofsLocal{subId};
            
            % SENSE SORT: Agrupem mantenint l'ordre estricte
            remLoc = [rLoc; dLoc];
            
            Krr = kMat(remLoc, remLoc);
            Krp = kMat(remLoc, pLoc);
            Kpr = kMat(pLoc, remLoc); 
            Kpp = kMat(pLoc, pLoc);
            
            fR = fVec(remLoc);
            fP = fVec(pLoc);
            
            numD = length(dLoc); 
            numR = length(remLoc);      
            
            % Com que hem concatenat [rLoc; dLoc], els duals estan al final
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
            coords = obj.meshCoords;
            gMinX = min(coords(:,1)); gMaxX = max(coords(:,1));
            gMinY = min(coords(:,2)); gMaxY = max(coords(:,2));
            
            numP = length(allPrimals);
            isDirichlet = false(numP, 1);
            
            for k = 1:numP
                node = allPrimals(k);
                x = coords(node, 1);
                y = coords(node, 2);
                
                if (abs(x - gMinX) < obj.nodeTol) || (abs(x - gMaxX) < obj.nodeTol) || ...
                   (abs(y - gMinY) < obj.nodeTol) || (abs(y - gMaxY) < obj.nodeTol)
                    isDirichlet(k) = true;
                end
            end
            
            activeIdx = find(~isDirichlet);  
        end
    end
end