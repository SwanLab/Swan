classdef FetiDPPoissonTop < handle
    % FETIDPPOISSON Class for the algebraic logic of the FETI-DP method.
    % 100% Matrix-Free Implementation for extreme scalability.
    
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
        
        dualIdxLocal       
        dualSignsLocal     
        dualWeights        
        
        % --- OPTIMIZED MATRIX-FREE CACHE ---
        localKiiFact       
        localKid           
        localKdi           
        localKdd           
        
        % --- CACHE FOR F(lambda) MATRIX-FREE ---
        localFDualBlock    
        sppActive          
        brActive           
        
        % --- NEW: CACHE FOR INSTANT RECONSTRUCTION ---
        localUrp
        localUrd
        localKrrInvFr
        localBasePrimalRhs
        localTdrUrp
    end
    
    methods (Access = public)
        
        function obj = FetiDPPoissonTop(globalMesh, subMeshes, stiffness, forces, tol, dirDofs)
            obj.meshCoords     = globalMesh.coord;
            obj.localMeshes    = subMeshes;
            obj.localStiffness = stiffness;
            obj.localForces    = forces;
            obj.nodeTol        = tol;
            
            subdomainsVec      = [size(subMeshes, 2) size(subMeshes, 1)];
            obj.numSubdomains  = prod(subdomainsVec);
            obj.dirichletDofs  = dirDofs;
            
            obj.extractFetiDofs();
        end
        
        function dBar = assembleProblem(obj)
            allPrimals = unique(vertcat(obj.primalDofsGlobal{:}));
            allDuals   = setdiff(unique(vertcat(obj.dualDofsGlobal{:})), allPrimals);
            
            numPrimals = length(allPrimals);
            numDuals   = length(allDuals);
            
            dBar        = zeros(numDuals, 1);
            SPP         = sparse(numPrimals, numPrimals);
            BrKrrInvKrp = sparse(numDuals, numPrimals);
            rhsPrimal   = zeros(numPrimals, 1);
            
            visitedDuals = zeros(max(allDuals), 1); 
            
            obj.localFDualBlock    = cell(obj.numSubdomains, 1);
            obj.localUrp           = cell(obj.numSubdomains, 1);
            obj.localUrd           = cell(obj.numSubdomains, 1);
            obj.localKrrInvFr      = cell(obj.numSubdomains, 1);
            obj.localBasePrimalRhs = cell(obj.numSubdomains, 1);
            obj.localTdrUrp        = cell(obj.numSubdomains, 1);
                        
            for subId = 1:obj.numSubdomains
                [Krr, Krp, Kpr, Kpp, fR, fP, Tdr] = obj.splitLocalMatrices(subId);
                
                Urd = Krr \ full(Tdr');     
                Urp = Krr \ full(Krp);             
                Ap  = obj.createAp(subId, allPrimals);
                
                % Guardamos TODO en caché para que la reconstrucción sea instantánea
                obj.localUrd{subId}           = Urd;
                obj.localUrp{subId}           = Urp;
                obj.localKrrInvFr{subId}      = Krr \ fR;
                obj.localBasePrimalRhs{subId} = fP - Urp' * fR;
                obj.localTdrUrp{subId}        = full(Tdr * Urp);
                
                dGlobal       = obj.dualDofsGlobal{subId};
                [~, dualRows] = ismember(dGlobal, allDuals);
                isFirst       = ~visitedDuals(dGlobal);
                
                dualSigns           = ones(length(dGlobal), 1);
                dualSigns(~isFirst) = -1;
                visitedDuals(dGlobal) = true;
                
                obj.dualIdxLocal{subId}   = dualRows;
                obj.dualSignsLocal{subId} = dualSigns;
                                
                sppLocal = Kpp - Kpr * Urp;
                SPP      = SPP + Ap * sppLocal * Ap';
                
                if ~isempty(dGlobal)
                    obj.localFDualBlock{subId} = Urd' * full(Tdr');
                    
                    dBar(dualRows)           = dBar(dualRows) + dualSigns .* (Urd' * fR);
                    BrKrrInvKrp(dualRows, :) = BrKrrInvKrp(dualRows, :) + dualSigns .* obj.localTdrUrp{subId} * Ap';
                end
                
                rhsPrimal = rhsPrimal + Ap * obj.localBasePrimalRhs{subId};
            end
                        
            activeDofs = obj.getActivePrimalDofs(allPrimals);   
            
            obj.sppActive = SPP(activeDofs, activeDofs);                    
            obj.brActive  = BrKrrInvKrp(:, activeDofs);               
            rhsActive     = rhsPrimal(activeDofs);                     
            
            dBar = dBar - obj.brActive * (obj.sppActive \ rhsActive);
            
            obj.setupPreconditioner();
        end

        function fLambda = applyFMat(obj, lambda)
            numDuals = length(lambda);
            fLambda  = zeros(numDuals, 1);
            
            for subId = 1:obj.numSubdomains
                dualRows = obj.dualIdxLocal{subId};
                if isempty(dualRows), continue; end
                
                dualSigns = obj.dualSignsLocal{subId};
                lamLocal  = dualSigns .* lambda(dualRows);
                
                resLocal = obj.localFDualBlock{subId} * lamLocal;
                fLambda(dualRows) = fLambda(dualRows) + dualSigns .* resLocal;
            end
            
            primalTerm = obj.brActive * (obj.sppActive \ full(obj.brActive' * lambda));
            fLambda    = fLambda + primalTerm;
        end
        
        function uFull = reconstructGlobalSolution(obj, lambdaSol, numNodes)
            allPrimals = unique(vertcat(obj.primalDofsGlobal{:}));
            numPrimals = length(allPrimals);
            rhsPrimal  = zeros(numPrimals, 1);
            
            % 1. Extraer RHS Primal sin resolver nada nuevo
            for subId = 1:obj.numSubdomains
                Ap = obj.createAp(subId, allPrimals);
                
                dualRows    = obj.dualIdxLocal{subId};
                dualSigns   = obj.dualSignsLocal{subId};
                lambdaLocal = dualSigns .* lambdaSol(dualRows);
        
                term      = obj.localTdrUrp{subId}' * lambdaLocal;
                rhsPrimal = rhsPrimal + Ap * (obj.localBasePrimalRhs{subId} + term);
            end
            
            activeDofs = obj.getActivePrimalDofs(allPrimals);
            uP         = zeros(numPrimals, 1);
            uP(activeDofs) = obj.sppActive \ rhsPrimal(activeDofs);
            
            % 2. Reconstrucción hiper-rápida (Zero matrices inversas)
            uFull = zeros(numNodes, 1); 
            
            for subId = 1:obj.numSubdomains
                pGlobal = obj.primalDofsGlobal{subId};
                dGlobal = obj.dualDofsGlobal{subId};
                rGlobal = obj.remDofsGlobal{subId};
                
                Ap          = obj.createAp(subId, allPrimals);
                dualRows    = obj.dualIdxLocal{subId};
                dualSigns   = obj.dualSignsLocal{subId};
                lambdaLocal = dualSigns .* lambdaSol(dualRows);
      
                uPLocal   = Ap' * uP;
                uRemLocal = obj.localKrrInvFr{subId} - obj.localUrp{subId} * uPLocal - obj.localUrd{subId} * lambdaLocal;
                
                uFull(pGlobal)      = uPLocal;
                allRemGlobal        = [rGlobal; dGlobal];
                uFull(allRemGlobal) = uRemLocal;
            end
        end
        
        function z = applyDirichletPrecond(obj, r)
            numDuals = length(r);
            z        = zeros(numDuals, 1);
            
            for subId = 1:obj.numSubdomains
                dualRows  = obj.dualIdxLocal{subId};
                if isempty(dualRows), continue; end
                
                dualSigns = obj.dualSignsLocal{subId};
                w         = obj.dualWeights(dualRows);   
                rLocal    = dualSigns .* r(dualRows);
                
                kidLocal = obj.localKid{subId};
                kdiLocal = obj.localKdi{subId};
                kddLocal = obj.localKdd{subId};
                kiiFact  = obj.localKiiFact{subId};
                
                z1 = kidLocal * (w .* rLocal);
                z2 = kiiFact \ z1; 
                z3 = kdiLocal * z2;
                z4 = kddLocal * (w .* rLocal) - z3;
                
                z(dualRows) = z(dualRows) + dualSigns .* (w .* z4);
            end
        end
        
        function n = getNumDuals(obj)
            allPrimals = unique(vertcat(obj.primalDofsGlobal{:}));
            allDuals   = setdiff(unique(vertcat(obj.dualDofsGlobal{:})), allPrimals);
            n          = length(allDuals);
        end
    end
    
    methods (Access = private)
        
        function extractFetiDofs(obj)
            nSub = prod(obj.numSubdomains);
            tol  = obj.nodeTol;
            
            pGlobal = cell(nSub,1); dGlobal = cell(nSub,1); rGlobal = cell(nSub,1);
            pLocal  = cell(nSub,1); dLocal  = cell(nSub,1); rLocal  = cell(nSub,1);
            
            for i = 1:nSub
                coords = obj.localMeshes{i}.coord;
                connec = obj.localMeshes{i}.connec;
                nNodes = size(coords, 1);
                
                [~, globalNodes] = ismembertol(coords, obj.meshCoords, tol, 'ByRows', true);
                
                allEdges       = [connec(:,[1,2]); connec(:,[2,3]); connec(:,[3,1])];
                allEdgesSorted = sort(allEdges, 2);
                [~, ia, ic]    = unique(allEdgesSorted, 'rows');
                edgeCount      = accumarray(ic, 1);
                boundaryEdges  = allEdges(ia(edgeCount == 1), :);
                boundaryLocalIdx = unique(boundaryEdges(:));
                
                bMeshes = obj.localMeshes{i}.createBoundaryMesh();
                
                primalGlobalSet = [];
                for b = 1:length(bMeshes)
                    gc          = bMeshes{b}.globalConnec;  
                    allN        = gc(:);
                    uN          = unique(allN);
                    deg         = accumarray(allN, 1, [max(allN), 1]);
                    extremLocal = uN(deg(uN) == 1);
                    primalGlobalSet = [primalGlobalSet; globalNodes(extremLocal)];
                end
                primalGlobalSet = unique(primalGlobalSet);
                
                [~, primalLocalIdx] = ismember(primalGlobalSet, globalNodes);
                primalLocalIdx      = primalLocalIdx(primalLocalIdx > 0);
                
                bCoords = coords(boundaryLocalIdx, :);
                minX = min(bCoords(:,1));  maxX = max(bCoords(:,1));
                minY = min(bCoords(:,2));  maxY = max(bCoords(:,2));
                
                isOnExtreme = (abs(coords(:,1) - minX) < tol) | ...
                              (abs(coords(:,1) - maxX) < tol) | ...
                              (abs(coords(:,2) - minY) < tol) | ...
                              (abs(coords(:,2) - maxY) < tol);
                              
                extremeLocalIdx = find(isOnExtreme);
                dualLocalIdx    = setdiff(intersect(extremeLocalIdx, boundaryLocalIdx), primalLocalIdx);
                
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
            obj.remDofsGlobal    = rGlobal; obj.remDofsLocal    = rLocal;
        end
        
        function [Krr, Krp, Kpr, Kpp, fR, fP, Tdr] = splitLocalMatrices(obj, subId)
            kMat = obj.localStiffness{subId};
            fVec = obj.localForces{subId};
            
            pLoc   = obj.primalDofsLocal{subId};
            dLoc   = obj.dualDofsLocal{subId};
            rLoc   = obj.remDofsLocal{subId};
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
            cols      = (1:numPrimalsLocal)';
            vals      = ones(numPrimalsLocal, 1);
            
            ap = sparse(rows, cols, vals, numPrimalsGlobal, numPrimalsLocal);
        end
        
        function activeIdx = getActivePrimalDofs(obj, allPrimals)
            isDirichlet = ismember(allPrimals, obj.dirichletDofs);
            activeIdx   = find(~isDirichlet);  
        end
        
        function [Kii, Kid, Kdi, Kdd] = computeLocalMatrices(obj, subId)
            kMat = obj.localStiffness{subId};
            rLoc = obj.remDofsLocal{subId};   
            dLoc = obj.dualDofsLocal{subId};  
        
            Kii = kMat(rLoc, rLoc);
            Kid = kMat(rLoc, dLoc);
            Kdi = kMat(dLoc, rLoc);
            Kdd = kMat(dLoc, dLoc);
        end
        
        function setupPreconditioner(obj)
            numDuals     = obj.getNumDuals();
            multiplicity = zeros(numDuals, 1);
            
            obj.localKiiFact     = cell(obj.numSubdomains, 1);
            obj.localKid         = cell(obj.numSubdomains, 1);
            obj.localKdi         = cell(obj.numSubdomains, 1);
            obj.localKdd         = cell(obj.numSubdomains, 1);
            
            for subId = 1:obj.numSubdomains
                [Kii, Kid, Kdi, Kdd] = obj.computeLocalMatrices(subId);
                
                obj.localKid{subId} = Kid;
                obj.localKdi{subId} = Kdi;
                obj.localKdd{subId} = Kdd;
                
                obj.localKiiFact{subId} = decomposition(Kii, 'chol');
                
                dualRows               = obj.dualIdxLocal{subId};
                multiplicity(dualRows) = multiplicity(dualRows) + 1;
            end
            
            obj.dualWeights = 1 ./ multiplicity;
        end
    end
end