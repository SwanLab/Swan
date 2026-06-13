classdef FetiDPPoissonV2 < handle
    % FETIDPPOISSON Classe per a la lògica algebraica de FETI-DP.
    
    properties (Access = private)
        localStiffness  
        localForces     
        localMeshes     
        meshCoords      
        nodeTol         
        numSubdomains   
        dirichletDofs  
        
        primalDofsGlobal      
        dualDofsGlobal        
        remDofsGlobal 
        
        primalDofsLocal    
        dualDofsLocal        
        remDofsLocal
        dualIdxLocal   % cell: índexs globals dins allDuals per cada subdomini
        dualSignsLocal % cell: signes (+1/-1)

        localSchurBlocks  % Cell array Sdd
        dualWeights
    end
    
    properties (Access = public)         
    end
    
    methods (Access = public)
        function obj = FetiDPPoissonV2(globalMesh, subMeshes, stiffness, forces, tol, dirDofs)
            obj.meshCoords     = globalMesh.coord;
            obj.localMeshes    = subMeshes;
            obj.localStiffness = stiffness;
            obj.localForces    = forces;
            obj.nodeTol        = tol;
            subdomainsVec      = [size(subMeshes, 2) size(subMeshes, 1)];
            obj.numSubdomains  = prod(subdomainsVec);
            obj.dirichletDofs  = dirDofs;
            %obj.bdLocal        = cell(obj.numSubdomains, 1);
            
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
                
                Ap = obj.createAp(subId, allPrimals);
                % [Bd, visitedDuals] = obj.createBd(subId, allDuals, visitedDuals);
                % obj.bdLocal{subId} = Bd;
                % Alternativa a Bd
                dGlobal = obj.dualDofsGlobal{subId};
                [~, dualRows] = ismember(dGlobal, allDuals);
                isFirst    = ~visitedDuals(dGlobal);
                dualSigns  = ones(length(dGlobal), 1);
                dualSigns(~isFirst) = -1;
                visitedDuals(dGlobal) = true;
                % Guardem per la reconstrucció (substitueix bdLocal)
                obj.dualIdxLocal{subId}   = dualRows;
                obj.dualSignsLocal{subId} = dualSigns;
                                
                sppLoc = Kpp - Kpr * Urp;
                SPP = SPP + Ap * sppLoc * Ap';
                
                if ~isempty(dGlobal)
                    %fDual = fDual + Bd * (Urd' * Tdr') * Bd';
                    fDual(dualRows, dualRows) = fDual(dualRows, dualRows) + (dualSigns) .* (Urd'*Tdr') .* dualSigns';
                    
                    %dBar  = dBar + Bd * (Urd' * fR);
                    dBar(dualRows) = dBar(dualRows) + dualSigns .* (Urd' * fR);
                    % BrKrrInvKrp = BrKrrInvKrp + Bd * (Tdr * Urp) * Ap';
                    BrKrrInvKrp(dualRows, :) = BrKrrInvKrp(dualRows, :) + dualSigns .* (Tdr * Urp) * Ap';
                end
                
                rhsPrimal = rhsPrimal + Ap * (fP - Urp' * fR);
            end
                        
            activeDofs = obj.getActivePrimalDofs(allPrimals);   
            
            SPPActive = SPP(activeDofs, activeDofs);                    
            BrActive  = BrKrrInvKrp(:, activeDofs);               
            rhsActive = rhsPrimal(activeDofs);                     
            
            fMat = fDual + BrActive * (SPPActive \ full(BrActive'));    
            dBar = dBar  - BrActive * (SPPActive \ rhsActive);

            obj.setupPreconditioner();
        end
        
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
                Ap = obj.createAp(subId, allPrimals);
                % Bd = obj.bdLocal{subId};
                
                Urp = Krr \ full(Krp); 
                Urd = Krr \ full(Tdr');
                KrrInvFr = Krr \ fR;
                
                UrpCell{subId} = Urp;
                UrdCell{subId} = Urd;
                KrrInvFrCell{subId} = KrrInvFr;
                
                sppLoc = Kpp - Kpr * Urp;
                SPP = SPP + Ap * sppLoc * Ap';
                
                %term = Tdr' * (Bd' * lambdaSol) - fR;
                dualRows  = obj.dualIdxLocal{subId};
                dualSigns = obj.dualSignsLocal{subId};
                lambdaLoc = dualSigns .* lambdaSol(dualRows);
        
                term = Tdr' * lambdaLoc - fR;
                rhsPrimal = rhsPrimal + Ap * (fP + Urp' * term);
            end
            
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
                %Bd = obj.bdLocal{subId};
                dualRows  = obj.dualIdxLocal{subId};
                dualSigns = obj.dualSignsLocal{subId};
                lambdaLoc = dualSigns .* lambdaSol(dualRows);
      
                uPLoc = Ap' * uP;
                %uRemLoc = KrrInvFr - Urp * uPLoc - Urd * (Bd' * lambdaSol);
                uRemLoc  = KrrInvFrCell{subId} - UrpCell{subId} * uPLoc ...
                   - UrdCell{subId} * lambdaLoc;
                
                uFull(pGlobal) = uPLoc;
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
            
            figure('Name', 'FETI-DP Nodes (Poisson)', 'Color', 'w');
            hold on; axis equal;
            for i = 1:obj.numSubdomains
                locCoords = obj.localMeshes{i}.coord;
                locConnec = obj.localMeshes{i}.connec;                
                patch('Faces', locConnec, 'Vertices', locCoords, ...
                      'FaceColor', 'none', 'EdgeColor', [0.8 0.8 0.8], ...
                      'LineWidth', 1.5, 'HandleVisibility', 'off');
            end
            
            scatter(rCoords(:,1), rCoords(:,2), 20, [0.5 0.5 0.5], 'filled', 'DisplayName', 'Internal (Remaining)');
            scatter(dCoords(:,1), dCoords(:,2), 40, 'b', 'filled', 'DisplayName', 'Interface + Boundary (Dual)');
            scatter(pCoords(:,1), pCoords(:,2), 80, 'r', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1, 'DisplayName', 'Corners (Primal)');
            
            legend('Location', 'bestoutside');
            title('FETI-DP Nodes Distribution (Poisson)');
            xlabel('X'); ylabel('Y');
            grid on;
            hold off;
        end

        % function [z, M] = applyDirichletPrecond(obj, r, fMat)
        %     allPrimals = unique(vertcat(obj.primalDofsGlobal{:}));
        %     allDuals   = setdiff(unique(vertcat(obj.dualDofsGlobal{:})), allPrimals);
        %     numDuals   = length(allDuals);
        % 
        %     % Calcula multiplicitat de cada dual global
        %     multiplicity = zeros(numDuals, 1);
        %     for subId = 1:obj.numSubdomains
        %         dualRows = obj.dualIdxLocal{subId};
        %         multiplicity(dualRows) = multiplicity(dualRows) + 1;
        %     end
        % 
        %     z = zeros(numDuals, 1);
        %     M = zeros(numDuals);
        % 
        %     for subId = 1:obj.numSubdomains
        %         Sdd = obj.computeLocalSchur(subId);   
        % 
        %         dualRows  = obj.dualIdxLocal{subId};
        %         dualSigns = obj.dualSignsLocal{subId};
        % 
        %         if isempty(dualRows), continue; end
        % 
        %         % Pesos locals
        %         w = 1./multiplicity(dualRows);
        % 
        %         % Extreu la part local de r (amb signe de Bd^T)
        %         rLoc = dualSigns .* r(dualRows);
        % 
        %         % Aplica el Schur local i escampa (amb signe de Bd)
        %         z(dualRows) = z(dualRows) + dualSigns .* (w .* (Sdd *(w .* rLoc)));
        % 
        %         D = diag(w);
        %         M(dualRows,dualRows) = M(dualRows,dualRows) + D * Sdd * D;
        %     end
        % end

        % Mètode per aplicar el precondicionador (interfície per PCG, 1 sol argument)
        % function z = applyDirichletPrecond(obj, r)
        %     [z, ~] = obj.applyDirichletPrecondFull(r);
        % end
        
        % % Mètode complet que retorna també la matriu M (per calcular κ)
        % function [z, M] = applyDirichletPrecondFull(obj, r)
        %     allPrimals = unique(vertcat(obj.primalDofsGlobal{:}));
        %     allDuals   = setdiff(unique(vertcat(obj.dualDofsGlobal{:})), allPrimals);
        %     numDuals   = length(allDuals);
        % 
        %     multiplicity = zeros(numDuals, 1);
        %     for subId = 1:obj.numSubdomains
        %         dualRows = obj.dualIdxLocal{subId};
        %         multiplicity(dualRows) = multiplicity(dualRows) + 1;
        %     end
        % 
        %     z = zeros(numDuals, 1);
        %     M = zeros(numDuals);
        % 
        %     for subId = 1:obj.numSubdomains
        %         Sdd       = obj.computeLocalSchur(subId);
        %         dualRows  = obj.dualIdxLocal{subId};
        %         dualSigns = obj.dualSignsLocal{subId};
        %         if isempty(dualRows), continue; end
        % 
        %         w    = 1 ./ multiplicity(dualRows);
        %         rLoc = dualSigns .* r(dualRows);
        %         z(dualRows) = z(dualRows) + dualSigns .* (w .* (Sdd * (w .* rLoc)));
        % 
        %         D = diag(w);
        %         M(dualRows, dualRows) = M(dualRows, dualRows) + D * Sdd * D;
        %     end
        % end

        % function [z, M] = applyDirichletPrecondFull(obj, r)
        %     numDuals = length(r);
        %     z = zeros(numDuals, 1);
        % 
        %     M = zeros(numDuals); 
        % 
        %     for subId = 1:obj.numSubdomains
        %         dualRows  = obj.dualIdxLocal{subId};
        %         dualSigns = obj.dualSignsLocal{subId};
        %         if isempty(dualRows), continue; end
        % 
        %         Sdd  = obj.localSchurBlocks{subId}; 
        %         w    = obj.dualWeights(dualRows);   
        % 
        %         rLoc = dualSigns .* r(dualRows);
        %         z(dualRows) = z(dualRows) + dualSigns .* (w .* (Sdd * (w .* rLoc)));
        % 
        %         % Opcional: Solo si de verdad necesitas construir M
        %         D = diag(w);
        %         M(dualRows, dualRows) = M(dualRows, dualRows) + D * Sdd * D;
        %     end
        % end
        % 
        % % Mètode públic per obtenir M sense aplicar-lo a cap vector
        % function M = buildPrecondMatrix(obj)
        %     [~, M] = obj.applyDirichletPrecondFull(zeros(obj.getNumDuals(), 1));
        % end

        % =========================================================
        % PRECONDICIONADOR DE DIRICHLET
        % =========================================================
        
        function z = applyDirichletPrecond(obj, r)
            % Aplica el precondicionador al vector r iterativamente (súper rápido).
            % Esta función es la que se debe pasar al solver PCG.
            numDuals = length(r);
            z = zeros(numDuals, 1);
            
            for subId = 1:obj.numSubdomains
                dualRows  = obj.dualIdxLocal{subId};
                if isempty(dualRows), continue; end
                
                dualSigns = obj.dualSignsLocal{subId};
                Sdd       = obj.localSchurBlocks{subId}; 
                w         = obj.dualWeights(dualRows);   
                
                % Multiplicación vectorizada (equivalente a Bd * D * Sdd * D * Bd^T * r)
                rLoc = dualSigns .* r(dualRows);
                z(dualRows) = z(dualRows) + dualSigns .* (w .* (Sdd * (w .* rLoc)));
            end
        end
        
        function M = buildPrecondMatrix(obj)
            % Construye la matriz explícita M del precondicionador.
            % Usar SOLO para análisis (ej. calcular condest(M * fMat)).
            numDuals = obj.getNumDuals();
            M = zeros(numDuals); 
            
            for subId = 1:obj.numSubdomains
                dualRows  = obj.dualIdxLocal{subId};
                if isempty(dualRows), continue; end
                
                dualSigns = obj.dualSignsLocal{subId};
                Sdd       = obj.localSchurBlocks{subId}; 
                w         = obj.dualWeights(dualRows);   
                
                % Incluimos el signo (Bd) y los pesos (w) en una matriz diagonal
                W_signed = diag(w .* dualSigns);
                
                % Ensamblaje aditivo de la matriz M
                M(dualRows, dualRows) = M(dualRows, dualRows) + W_signed * Sdd * W_signed;
            end
        end
               
        function n = getNumDuals(obj)
            allPrimals = unique(vertcat(obj.primalDofsGlobal{:}));
            allDuals   = setdiff(unique(vertcat(obj.dualDofsGlobal{:})), allPrimals);
            n = length(allDuals);
        end
    end
    
    methods (Access = private)
        % function extractFetiDofs(obj)
        %     nSub = prod(obj.numSubdomains);
        %     tol = obj.nodeTol;
        % 
        %     pGlobal = cell(nSub, 1); dGlobal = cell(nSub, 1); rGlobal = cell(nSub, 1);
        %     pLocal = cell(nSub, 1);  dLocal = cell(nSub, 1);  rLocal = cell(nSub, 1);
        % 
        %     for i = 1:nSub
        %         coords = obj.localMeshes{i}.coord;
        % 
        %         minX = min(coords(:,1)); maxX = max(coords(:,1));
        %         minY = min(coords(:,2)); maxY = max(coords(:,2));
        % 
        %         isBL = abs(coords(:,1) - minX) < tol & abs(coords(:,2) - minY) < tol;
        %         isBR = abs(coords(:,1) - maxX) < tol & abs(coords(:,2) - minY) < tol;
        %         isTL = abs(coords(:,1) - minX) < tol & abs(coords(:,2) - maxY) < tol;
        %         isTR = abs(coords(:,1) - maxX) < tol & abs(coords(:,2) - maxY) < tol;
        % 
        %         localPrimalNodes = find(isBL | isBR | isTL | isTR);
        % 
        %         isOnLocalBoundary = abs(coords(:,1) - minX) < tol | abs(coords(:,1) - maxX) < tol | ...
        %                             abs(coords(:,2) - minY) < tol | abs(coords(:,2) - maxY) < tol;
        % 
        %         isDual = isOnLocalBoundary; 
        %         localDualNodes = setdiff(find(isDual), localPrimalNodes);
        % 
        %         isRem = ~isDual;
        %         localRemNodes = setdiff(find(isRem), localPrimalNodes);
        % 
        %         [~, globalNodes] = ismembertol(coords, obj.meshCoords, tol, 'ByRows', true);                
        % 
        %         pGlobal{i} = globalNodes(localPrimalNodes);
        %         dGlobal{i} = globalNodes(localDualNodes);
        %         rGlobal{i} = globalNodes(localRemNodes);
        % 
        %         pLocal{i} = localPrimalNodes;
        %         dLocal{i} = localDualNodes;
        %         rLocal{i} = localRemNodes;
        %     end
        %     obj.primalDofsLocal = pLocal;
        %     obj.dualDofsLocal   = dLocal;
        %     obj.remDofsLocal    = rLocal;
        % 
        %     obj.primalDofsGlobal = pGlobal; 
        %     obj.dualDofsGlobal   = dGlobal; 
        %     obj.remDofsGlobal    = rGlobal;    
        % end

        function extractFetiDofs(obj)
            nSub = prod(obj.numSubdomains);
            tol  = obj.nodeTol;

            pGlobal = cell(nSub,1); 
            dGlobal = cell(nSub,1); 
            rGlobal = cell(nSub,1);
            pLocal  = cell(nSub,1); 
            dLocal  = cell(nSub,1); 
            rLocal  = cell(nSub,1);

            for i = 1:nSub
                coords = obj.localMeshes{i}.coord;
                connec = obj.localMeshes{i}.connec;
                nNodes = size(coords, 1);

                % Mapeig local a global
                [~, globalNodes] = ismembertol(coords, obj.meshCoords, tol, 'ByRows', true);

                % Extreure arestes de contorn del subdomini
                allEdges       = [connec(:,[1,2]); connec(:,[2,3]); connec(:,[3,1])];
                allEdgesSorted = sort(allEdges, 2);
                [~, ia, ic]    = unique(allEdgesSorted, 'rows');
                edgeCount      = accumarray(ic, 1);
                boundaryEdges  = allEdges(ia(edgeCount == 1), :);
                boundaryLocalIdx = unique(boundaryEdges(:));

               % Obtenir les boundary meshes del subdomini
                bMeshes = obj.localMeshes{i}.createBoundaryMesh();

                % PRIMAL: extrems de cada boundary mesh
                primalGlobalSet = [];
                for b = 1:length(bMeshes)
                    gc = bMeshes{b}.globalConnec;  

                    allN = gc(:);
                    uN   = unique(allN);
                    deg  = accumarray(allN, 1, [max(allN), 1]);

                    extremLocal = uN(deg(uN) == 1);
                    primalGlobalSet = [primalGlobalSet; globalNodes(extremLocal)];
                end
                primalGlobalSet = unique(primalGlobalSet);

                % Converteix a índexs locals
                [~, primalLocalIdx] = ismember(primalGlobalSet, globalNodes);
                primalLocalIdx = primalLocalIdx(primalLocalIdx > 0);

               % DUAL: nodes en les 4 línies extremes del bounding box - primaris
                bCoords = coords(boundaryLocalIdx, :);
                minX = min(bCoords(:,1));  maxX = max(bCoords(:,1));
                minY = min(bCoords(:,2));  maxY = max(bCoords(:,2));

                isOnExtreme = (abs(coords(:,1) - minX) < tol) | ...
                              (abs(coords(:,1) - maxX) < tol) | ...
                              (abs(coords(:,2) - minY) < tol) | ...
                              (abs(coords(:,2) - maxY) < tol);

                extremeLocalIdx = find(isOnExtreme);
                dualLocalIdx    = setdiff(intersect(extremeLocalIdx, boundaryLocalIdx), primalLocalIdx);

                % REMAINING
                remLocalIdx = setdiff((1:nNodes)', [primalLocalIdx; dualLocalIdx]);

                pGlobal{i} = globalNodes(primalLocalIdx);
                dGlobal{i} = globalNodes(dualLocalIdx);
                rGlobal{i} = globalNodes(remLocalIdx);

                pLocal{i}  = primalLocalIdx;
                dLocal{i}  = dualLocalIdx;
                rLocal{i}  = remLocalIdx;
            end

            obj.primalDofsGlobal = pGlobal; obj.primalDofsLocal = pLocal;
            obj.dualDofsGlobal   = dGlobal; obj.dualDofsLocal   = dLocal;
            obj.remDofsGlobal    = rGlobal; obj.remDofsLocal     = rLocal;
        end
        
        function [Krr, Krp, Kpr, Kpp, fR, fP, Tdr] = splitLocalMatrices(obj, subId)
            kMat = obj.localStiffness{subId};
            fVec = obj.localForces{subId};
            
            pLoc = obj.primalDofsLocal{subId};
            dLoc = obj.dualDofsLocal{subId};
            rLoc = obj.remDofsLocal{subId};
            
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
            
            if numD == 0
                rows = []; cols = []; vals = [];
            else
                vals = ones(numD, 1);
            end
            
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
        
        % function [bd, visitedDuals] = createBd(obj, subId, allDuals, visitedDuals)
        %     dGlobal = obj.dualDofsGlobal{subId};
        %     numDualsGlobal = length(allDuals);
        %     numDualsLocal  = length(dGlobal);
        % 
        %     if numDualsLocal == 0
        %         bd = sparse(numDualsGlobal, 0);
        %         return;
        %     end
        % 
        %     [~, rows] = ismember(dGlobal, allDuals);
        %     cols = (1:numDualsLocal)';
        % 
        %     isFirstVisit = (visitedDuals(dGlobal) == 0);
        %     vals = ones(numDualsLocal, 1);
        %     vals(~isFirstVisit) = -1;
        % 
        %     visitedDuals(dGlobal(isFirstVisit)) = 1;
        % 
        %     bd = sparse(rows, cols, vals, numDualsGlobal, numDualsLocal);
        % end
           
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
        
            Sdd = Kdd - Kdi * (Kii \ Kid);
        end

        function setupPreconditioner(obj)
            numDuals = obj.getNumDuals();
            multiplicity = zeros(numDuals, 1);
            obj.localSchurBlocks = cell(obj.numSubdomains, 1);
            
            for subId = 1:obj.numSubdomains
                obj.localSchurBlocks{subId} = obj.computeLocalSchur(subId);
                
                dualRows = obj.dualIdxLocal{subId};
                multiplicity(dualRows) = multiplicity(dualRows) + 1;
            end

            obj.dualWeights = 1 ./ multiplicity;
        end
    end
end