classdef FetiDP < handle
    % FETIDP Classe encarregada de la lògica algebraica del mètode FETI-DP.
    % Calcula els graus de llibertat, ensambla el problema dual i reconstrueix la solució.
    
    properties (Access = private)
        localStiffness  
        localForces    
        localMeshes   
        meshCoords      
        dofsPerNode     % Graus de llibertat per node (1 per Poisson, 2 per Elasticitat)
        nodeTol         % Tolerància geomètrica per identificar nodes
        numSubdomains   
        dirichletDofs  
        
        primalDofs      % Graus de llibertat primals globals per subdomini
        dualDofs        % Graus de llibertat duals globals per subdomini
        remDofs         % Graus de llibertat restants globals per subdomini
    end
    
    properties (Access = public)
        bdLocal         % Cel·la que guarda les matrius Bd per subdomini
    end
    
    methods (Access = public)
        function obj = FetiDP(globalMesh, subMeshes, stiffness, forces, tol, dofsNode, dirDofs)
            % Constructor de la classe FetiDP
            obj.meshCoords     = globalMesh.coord;
            obj.dofsPerNode    = dofsNode;
            obj.localMeshes    = subMeshes;
            obj.localStiffness = stiffness;
            obj.localForces    = forces;
            obj.nodeTol        = tol;
            obj.numSubdomains  = length(stiffness);
            obj.dirichletDofs  = dirDofs;
            obj.bdLocal        = cell(obj.numSubdomains, 1);
            
            % Classificació dels graus de llibertat
            obj.extractFetiDofs();
        end
        
        function [fMat, dBar] = assembleProblem(obj)
            % ASSEMBLEPROBLEM Construeix la matriu del problema dual F i el vector dbar
            
            allPrimals = unique(vertcat(obj.primalDofs{:}));
            allDuals   = setdiff(unique(vertcat(obj.dualDofs{:})), allPrimals);
            
            numPrimals = length(allPrimals);
            numDuals   = length(allDuals);
            
            fDual       = sparse(numDuals, numDuals);
            dBar        = zeros(numDuals, 1);
            spp         = sparse(numPrimals, numPrimals);
            BrKrrInvKrp = sparse(numDuals, numPrimals);
            rhsPrimal   = zeros(numPrimals, 1);
            
            visitedDuals = zeros(max(allDuals), 1); 
            
            for subId = 1:obj.numSubdomains
                [kRr, kRp, kPr, kPp, fR, fP, tDr] = obj.splitLocalMatrices(subId);
                
                uRd = kRr \ full(tDr');     
                uRp = kRr \ full(kRp);             
                
                ap = obj.createAp(subId, allPrimals);
                [bd, visitedDuals] = obj.createBd(subId, allDuals, visitedDuals);
                obj.bdLocal{subId} = bd;
                                
                sppLoc = kPp - kPr * uRp;
                spp = spp + ap * sppLoc * ap';
                
                if size(bd, 2) > 0
                    fDual = fDual + bd * (uRd' * tDr') * bd';
                    dBar  = dBar + bd * (uRd' * fR);
                    BrKrrInvKrp = BrKrrInvKrp + bd * (tDr * uRp) * ap';
                end
                
                rhsPrimal = rhsPrimal + ap * (fP - uRp' * fR);
            end
            
            activeDofs = obj.getActivePrimalDofs(allPrimals);   
            
            sppActive = spp(activeDofs, activeDofs);                    
            brActive  = BrKrrInvKrp(:, activeDofs);               
            rhsActive = rhsPrimal(activeDofs);                     
            
            fMat = fDual + brActive * (sppActive \ full(brActive'));    
            dBar = dBar  - brActive * (sppActive \ rhsActive);         
        end
        
        function zVec = applyDirichletPrecond(obj, rVec)
            % APPLYDIRICHLETPRECOND Aplica el precondicionador de Dirichlet
            numL = length(rVec);
            zVec = zeros(numL, 1);
            scaling = 0.5; 
            
            for subId = 1:obj.numSubdomains
                bd = obj.bdLocal{subId};
                if size(bd, 2) == 0
                    continue; 
                end
                
                kMat = obj.localStiffness{subId};
                
                dGlobal = obj.dualDofs{subId};
                iGlobal = obj.remDofs{subId}; 
                
                allDofs = sort([obj.primalDofs{subId}; dGlobal; iGlobal]);
                [~, dLoc] = ismember(dGlobal, allDofs);
                [~, iLoc] = ismember(iGlobal, allDofs);
                
                kIi = kMat(iLoc, iLoc);
                kId = kMat(iLoc, dLoc);
                kDi = kMat(dLoc, iLoc);
                kDd = kMat(dLoc, dLoc);
                
                sLocal = kDd - kDi * (kIi \ kId);
                rLocal = scaling * (bd' * rVec);
                zLocal = sLocal * rLocal;
                
                zVec = zVec + bd * (scaling * zLocal);
            end
        end
        
        function uFull = reconstructGlobalSolution(obj, lambdaSol, numNodes)
            % RECONSTRUCTGLOBALSOLUTION Construeix el vector global de resultats
            
            allPrimals = unique(vertcat(obj.primalDofs{:}));
            numP = length(allPrimals);
            spp = sparse(numP, numP);
            rhsPrimal = zeros(numP, 1);
            
            UrpCell = cell(obj.numSubdomains, 1);
            UrdCell = cell(obj.numSubdomains, 1);
            KrrInvFrCell = cell(obj.numSubdomains, 1);
            
            for subId = 1:obj.numSubdomains
                [kRr, kRp, kPr, kPp, fR, fP, tDr] = obj.splitLocalMatrices(subId);
                ap = obj.createAp(subId, allPrimals);
                bd = obj.bdLocal{subId};
                
                uRp = kRr \ full(kRp); 
                uRd = kRr \ full(tDr');
                krrInvFr = kRr \ fR;
                
                UrpCell{subId} = uRp;
                UrdCell{subId} = uRd;
                KrrInvFrCell{subId} = krrInvFr;
                
                sppLoc = kPp - kPr * uRp;
                spp = spp + ap * sppLoc * ap';
                
                term = tDr' * (bd' * lambdaSol) - fR;
                rhsPrimal = rhsPrimal + ap * (fP + uRp' * term);
            end
            
            activeDofs = obj.getActivePrimalDofs(allPrimals);
            sppActive = spp(activeDofs, activeDofs);
            rhsActive = rhsPrimal(activeDofs);
            
            uP = zeros(numP, 1);
            uP(activeDofs) = sppActive \ rhsActive;
            
            uFull = zeros(numNodes * obj.dofsPerNode, 1);
            
            for subId = 1:obj.numSubdomains
                pGlobal = obj.primalDofs{subId};
                dGlobal = obj.dualDofs{subId};
                rGlobal = obj.remDofs{subId};
                
                ap = obj.createAp(subId, allPrimals);
                bd = obj.bdLocal{subId};
                
                uRp = UrpCell{subId};
                uRd = UrdCell{subId};
                krrInvFr = KrrInvFrCell{subId};
                
                uPLoc = ap' * uP;
                uRemLoc = krrInvFr - uRp * uPLoc - uRd * (bd' * lambdaSol);
                
                uFull(pGlobal) = uPLoc;
                
                allRemGlobalSorted = sort([rGlobal; dGlobal]);
                uFull(allRemGlobalSorted) = uRemLoc;
            end
        end
        
        function visualizeFetiNodes(obj)
            % VISUALIZEFETINODES Dibuixa l'esquema de nodes
            
            allPrimalDofs = unique(vertcat(obj.primalDofs{:}));
            allDualDofs   = unique(vertcat(obj.dualDofs{:}));
            allRemDofs    = unique(vertcat(obj.remDofs{:}));
            
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
    end
    
    methods (Access = private)
        function extractFetiDofs(obj)
            nSub = prod(obj.numSubdomains);
            pGlobal = cell(nSub, 1);
            dGlobal = cell(nSub, 1);
            rGlobal = cell(nSub, 1);

            tol = obj.nodeTol;
            ndim = obj.localMeshes{1,1}.ndim;
            globalCoords = obj.meshCoords; 

            % --- LIMITES GLOBALES DEL DOMINIO ---
            gMinX = min(globalCoords(:,1)); gMaxX = max(globalCoords(:,1));
            gMinY = min(globalCoords(:,2)); gMaxY = max(globalCoords(:,2));

            for i = 1:nSub
                coords = obj.localMeshes{i}.coord;

                % Limits locals del subdomini
                minX = min(coords(:,1)); maxX = max(coords(:,1));
                minY = min(coords(:,2)); maxY = max(coords(:,2));

                % Nodes Primals 
                isBL = abs(coords(:,1) - minX) < tol & abs(coords(:,2) - minY) < tol; % BottomLeft
                isBR = abs(coords(:,1) - maxX) < tol & abs(coords(:,2) - minY) < tol; % BottomRight
                isTL = abs(coords(:,1) - minX) < tol & abs(coords(:,2) - maxY) < tol; % TopLeft
                isTR = abs(coords(:,1) - maxX) < tol & abs(coords(:,2) - maxY) < tol; % TopRight

                localPrimalNodes = find(isBL | isBR | isTL | isTR);

                % Identificar tots els nodes a la frontera del subdomini
                isOnLocalBoundary = abs(coords(:,1) - minX) < tol | abs(coords(:,1) - maxX) < tol | ...
                                    abs(coords(:,2) - minY) < tol | abs(coords(:,2) - maxY) < tol;

                % Identificar fronteras globals lliures (Top, Bottom, Right)
                % Aquestes zones no han de ser duals per a que la malla
                % pugui deformarse allí
                isOnFreeGlobalBoundary = (abs(coords(:,1) - gMaxX) < tol) | ...
                                         (abs(coords(:,2) - gMinY) < tol) | ...
                                         (abs(coords(:,2) - gMaxY) < tol);

                % Nodes duals: Frontera local menys les fronteres lliures globals
                % Això deixa unicament les interficies inernes i la pared esquerra
                isDual = isOnLocalBoundary & ~isOnFreeGlobalBoundary;
                localDualNodes = setdiff(find(isDual), localPrimalNodes);

                % Nodes restants: Todo lo que no sea Dual (incluye interiores y fronteras libres)
                isRem = ~isDual;
                localRemNodes = setdiff(find(isRem), localPrimalNodes);

                % 6. Mapeo a índices globales y conversión a Grados de Libertad (DoFs)
                [~, globalNodes] = ismembertol(coords, globalCoords, tol, 'ByRows', true);                

                pGlobal{i} = obj.nodesToDofs(globalNodes(localPrimalNodes));
                dGlobal{i} = obj.nodesToDofs(globalNodes(localDualNodes));
                rGlobal{i} = obj.nodesToDofs(globalNodes(localRemNodes));
            end

            obj.primalDofs = pGlobal;
            obj.dualDofs   = dGlobal;
            obj.remDofs    = rGlobal;          
        end
        % function extractFetiDofs(obj)
        %     % EXTRACTFETIDOFS Classifica els graus de llibertat assegurant
        %     % que els límits globals NO es consideren duals.
        % 
        %     pGlobal = cell(obj.numSubdomains, 1);
        %     dGlobal = cell(obj.numSubdomains, 1);
        %     rGlobal = cell(obj.numSubdomains, 1);
        % 
        %     gMinX = min(obj.meshCoords(:,1)); gMaxX = max(obj.meshCoords(:,1));
        %     gMinY = min(obj.meshCoords(:,2)); gMaxY = max(obj.meshCoords(:,2));
        % 
        %     for i = 1:obj.numSubdomains
        %         coords = obj.localMeshes{i}.coord;
        % 
        %         minX = min(coords(:,1)); maxX = max(coords(:,1));
        %         minY = min(coords(:,2)); maxY = max(coords(:,2));
        % 
        %         isBL = abs(coords(:,1) - minX) < obj.nodeTol & abs(coords(:,2) - minY) < obj.nodeTol;
        %         isBR = abs(coords(:,1) - maxX) < obj.nodeTol & abs(coords(:,2) - minY) < obj.nodeTol;
        %         isTL = abs(coords(:,1) - minX) < obj.nodeTol & abs(coords(:,2) - maxY) < obj.nodeTol;
        %         isTR = abs(coords(:,1) - maxX) < obj.nodeTol & abs(coords(:,2) - maxY) < obj.nodeTol;
        % 
        %         localPrimalNodes = find(isBL | isBR | isTL | isTR);
        % 
        %         isOnLocalBoundary = abs(coords(:,1) - minX) < obj.nodeTol | abs(coords(:,1) - maxX) < obj.nodeTol | ...
        %                             abs(coords(:,2) - minY) < obj.nodeTol | abs(coords(:,2) - maxY) < obj.nodeTol;
        % 
        %         % Tota la frontera global (exterior) NO ha de ser Dual mai
        %         isOnGlobalBoundary = (abs(coords(:,1) - gMinX) < obj.nodeTol) | ...
        %                              (abs(coords(:,1) - gMaxX) < obj.nodeTol) | ...
        %                              (abs(coords(:,2) - gMinY) < obj.nodeTol) | ...
        %                              (abs(coords(:,2) - gMaxY) < obj.nodeTol);
        % 
        %         isDual = isOnLocalBoundary & ~isOnGlobalBoundary;
        %         localDualNodes = setdiff(find(isDual), localPrimalNodes);
        % 
        %         isRem = ~isDual;
        %         localRemNodes = setdiff(find(isRem), localPrimalNodes);
        % 
        %         [~, globalNodes] = ismembertol(coords, obj.meshCoords, obj.nodeTol, 'ByRows', true);                
        % 
        %         pGlobal{i} = obj.nodesToDofs(globalNodes(localPrimalNodes));
        %         dGlobal{i} = obj.nodesToDofs(globalNodes(localDualNodes));
        %         rGlobal{i} = obj.nodesToDofs(globalNodes(localRemNodes));
        %     end
        % 
        %     obj.primalDofs = pGlobal;
        %     obj.dualDofs   = dGlobal;
        %     obj.remDofs    = rGlobal;          
        % end
        
        function dofs = nodesToDofs(obj, nodes)            
            nodes = nodes(:);
            dofs = zeros(length(nodes) * obj.dofsPerNode, 1);
            for d = 1:obj.dofsPerNode
                dofs(d:obj.dofsPerNode:end) = (nodes-1) * obj.dofsPerNode + d;
            end
            dofs = sort(dofs);
        end
        
        function [kRr, kRp, kPr, kPp, fR, fP, tDr] = splitLocalMatrices(obj, subId)
            kMat = obj.localStiffness{subId};
            fVec = obj.localForces{subId};
            
            pDofsGlobal = obj.primalDofs{subId};
            dDofsGlobal = obj.dualDofs{subId};
            rDofsGlobal = obj.remDofs{subId};
            
            allDofsGlobal = sort([pDofsGlobal; dDofsGlobal; rDofsGlobal]);
            
            [~, pDofsLocal] = ismember(pDofsGlobal, allDofsGlobal);
            [~, dDofsLocal] = ismember(dDofsGlobal, allDofsGlobal);
            [~, rDofsLocal] = ismember(rDofsGlobal, allDofsGlobal);
            
            remDofsLocal = sort([rDofsLocal; dDofsLocal]);
            pDofsLocal   = sort(pDofsLocal);
            
            kRr = kMat(remDofsLocal, remDofsLocal);
            kRp = kMat(remDofsLocal, pDofsLocal);
            kPr = kMat(pDofsLocal, remDofsLocal); 
            kPp = kMat(pDofsLocal, pDofsLocal);
            
            fR = fVec(remDofsLocal);
            fP = fVec(pDofsLocal);
                        
            dDofsLocalSorted = sort(dDofsLocal);
            numD = length(dDofsLocalSorted); 
            numR = length(remDofsLocal);      
            
            [~, posInRem] = ismember(dDofsLocalSorted, remDofsLocal);
            
            rows = (1:numD)';
            cols = posInRem;
            vals = ones(numD, 1);
            
            tDr = sparse(rows, cols, vals, numD, numR);
        end
        
        function ap = createAp(obj, subId, allPrimals)
            pGlobal = obj.primalDofs{subId};
            numPrimalsGlobal = length(allPrimals);
            numPrimalsLocal  = length(pGlobal);
            
            [~, rows] = ismember(pGlobal, allPrimals);
            cols = (1:numPrimalsLocal)';
            vals = ones(numPrimalsLocal, 1);
            
            ap = sparse(rows, cols, vals, numPrimalsGlobal, numPrimalsLocal);
        end
        
        function [bd, visitedDuals] = createBd(obj, subId, allDuals, visitedDuals)
            dGlobal = obj.dualDofs{subId};
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
            % GETACTIVEPRIMALDOFS Filtra els primals utilitzant els dofs de
            % Dirichlet proporcionats per l'usuari a la inicialització.
            isDirichlet = ismember(allPrimals, obj.dirichletDofs);
            activeIdx = find(~isDirichlet);  
        end
    end
end