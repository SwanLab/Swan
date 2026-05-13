classdef FetiDPElasticity < handle
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
        
        % Preconditioner & Interface Assembly Data
        dualIdxLocal    
        dualSignsLocal  
    end
    
    % =========================================================
    % PUBLIC METHODS: INITIALIZATION & MAIN SOLVER
    % =========================================================
    methods (Access = public)
        
        function obj = FetiDPElasticity(globalMesh, subMeshes, stiffness, forces, tol, dofsNode, dirDofs)
            % Constructor: Initializes the FETI-DP solver and extracts DoFs
            obj.meshCoords     = globalMesh.coord;
            obj.dofsPerNode    = dofsNode;
            obj.localMeshes    = subMeshes;
            obj.localStiffness = stiffness;
            obj.localForces    = forces;
            obj.nodeTol        = tol;
            obj.numSubdomains  = length(stiffness);
            obj.dirichletDofs  = dirDofs;
            
            % Extract Primal, Dual, and Remaining DoFs
            obj.extractFetiDofs();
        end
        
        function [fMat, dBar] = assembleProblem(obj)
            % Assembles the global dual interface problem: fMat * lambda = dBar
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
                
                Urd = Krr \ full(Tdr');
                Urp = Krr \ full(Krp);
                Ap  = obj.createAp(subId, allPrimals);
                
                % Compute Boolean interface signs and local indices
                dGlobal       = obj.dualDofsGlobal{subId};
                [~, dualRows] = ismember(dGlobal, allDuals);
                
                isFirst             = ~visitedDuals(dGlobal);
                dualSigns           = ones(length(dGlobal), 1);
                dualSigns(~isFirst) = -1;
                visitedDuals(dGlobal) = true;
                
                % Store for reconstruction and preconditioner logic
                obj.dualIdxLocal{subId}   = dualRows;
                obj.dualSignsLocal{subId} = dualSigns;
                
                % Assemble Primal Schur Complement
                SppLoc = Kpp - Kpr * Urp;
                SPP    = SPP + Ap * SppLoc * Ap';
                
                % Assemble Dual Interface system
                if ~isempty(dGlobal)
                    fDual(dualRows, dualRows) = fDual(dualRows, dualRows) + dualSigns .* (Urd' * Tdr') .* dualSigns';
                    dBar(dualRows)            = dBar(dualRows) + dualSigns .* (Urd' * fR);
                    BrKrrInvKrp(dualRows, :)  = BrKrrInvKrp(dualRows, :) + dualSigns .* (Tdr * Urp) * Ap';
                end
                
                rhsPrimal = rhsPrimal + Ap * (fP - Urp' * fR);
            end
            
            % Reduce system to active (non-Dirichlet) Primal DoFs
            activeDofs = obj.getActivePrimalDofs(allPrimals);
            SPPActive  = SPP(activeDofs, activeDofs);
            BrActive   = BrKrrInvKrp(:, activeDofs);
            rhsActive  = rhsPrimal(activeDofs);
            
            % Final formulation of the dual interface matrix
            fMat = fDual + BrActive * (SPPActive \ BrActive');
            dBar = dBar  - BrActive * (SPPActive \ rhsActive);
        end
        
        function uFull = reconstructGlobalSolution(obj, lambdaSol, numNodes)
            % Reconstructs the full displacement field using the computed multipliers
            allPrimals = unique(vertcat(obj.primalDofsGlobal{:}));
            numP       = length(allPrimals);
            SPP        = sparse(numP, numP);
            rhsPrimal  = zeros(numP, 1);
            
            UrpCell      = cell(obj.numSubdomains, 1);
            UrdCell      = cell(obj.numSubdomains, 1);
            KrrInvFrCell = cell(obj.numSubdomains, 1);
            
            % Step 1: Accumulate primal contributions
            for subId = 1:obj.numSubdomains
                [Krr, Krp, Kpr, Kpp, fR, fP, Tdr] = obj.splitLocalMatrices(subId);
                Ap = obj.createAp(subId, allPrimals);
                
                Urp      = Krr \ Krp;
                Urd      = Krr \ Tdr';
                KrrInvFr = Krr \ fR;
                
                UrpCell{subId}      = Urp;
                UrdCell{subId}      = Urd;
                KrrInvFrCell{subId} = KrrInvFr;
                
                sppLoc    = Kpp - Kpr * Urp;
                SPP       = SPP + Ap * sppLoc * Ap';
                dualRows  = obj.dualIdxLocal{subId};
                dualSigns = obj.dualSignsLocal{subId};
                
                lambdaLoc = dualSigns .* lambdaSol(dualRows);
                term      = Tdr' * lambdaLoc - fR;
                rhsPrimal = rhsPrimal + Ap * (fP + Urp' * term);
            end
            
            % Solve for the active Primal DoFs
            activeDofs = obj.getActivePrimalDofs(allPrimals);
            sppActive  = SPP(activeDofs, activeDofs);
            rhsActive  = rhsPrimal(activeDofs);
            
            uP = zeros(numP, 1);
            uP(activeDofs) = sppActive \ rhsActive;
            
            % Step 2: Reconstruct local Remainder DoFs
            uFull = zeros(numNodes * obj.dofsPerNode, 1);
            
            for subId = 1:obj.numSubdomains
                pGlobal = obj.primalDofsGlobal{subId};
                dGlobal = obj.dualDofsGlobal{subId};
                rGlobal = obj.remDofsGlobal{subId};
                Ap      = obj.createAp(subId, allPrimals);
                
                dualRows  = obj.dualIdxLocal{subId};
                dualSigns = obj.dualSignsLocal{subId};
                lambdaLoc = dualSigns .* lambdaSol(dualRows);
                
                uPLoc   = Ap' * uP;
                uRemLoc = KrrInvFrCell{subId} - UrpCell{subId} * uPLoc - UrdCell{subId} * lambdaLoc;
                
                uFull(pGlobal)            = uPLoc;
                uFull([rGlobal; dGlobal]) = uRemLoc;
            end
        end
    end
    
    % =========================================================
    % PUBLIC METHODS: PRECONDITIONER
    % =========================================================
    methods (Access = public)
        
        function z = applyDirichletPrecond(obj, r)
            % Applies the Dirichlet preconditioner iteratively
            allPrimals = unique(vertcat(obj.primalDofsGlobal{:}));
            allDuals   = setdiff(unique(vertcat(obj.dualDofsGlobal{:})), allPrimals);
            numDuals   = length(allDuals);
            
            % Compute interface node multiplicity for scaling
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
            % Builds the explicit dense matrix M of the preconditioner for analysis
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
            % Visualizes the Domain Decomposition mapping (Primal, Dual, Internal)
            allPrimalDofs = unique(vertcat(obj.primalDofsGlobal{:}));
            allDualDofs   = unique(vertcat(obj.dualDofsGlobal{:}));
            allRemDofs    = unique(vertcat(obj.remDofsGlobal{:}));
            
            primalNodes = unique(ceil(allPrimalDofs / obj.dofsPerNode));
            dualNodes   = unique(ceil(allDualDofs   / obj.dofsPerNode));
            remNodes    = unique(ceil(allRemDofs    / obj.dofsPerNode));

            dualNodes = setdiff(dualNodes, primalNodes);
            remNodes  = setdiff(remNodes,  union(primalNodes, dualNodes));
            
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
            % Evaluates geometries and extracts Primal, Dual, and Internal DoFs
            nSub = obj.numSubdomains;
            tol  = obj.nodeTol;
            
            pGlobal = cell(nSub,1); dGlobal = cell(nSub,1); rGlobal = cell(nSub,1);
            pLocal  = cell(nSub,1); dLocal  = cell(nSub,1); rLocal  = cell(nSub,1);
            
            % Calculate node multiplicity to find shared interfaces
            nodeMultiplicity = zeros(size(obj.meshCoords,1), 1);
            for i = 1:nSub
                coords = obj.localMeshes{i}.coord;
                [~, gNodes] = ismembertol(coords, obj.meshCoords, tol, 'ByRows', true);
                nodeMultiplicity(gNodes) = nodeMultiplicity(gNodes) + 1;
            end
            
            % Map Dirichlet DoFs to Nodes
            isDirNode = false(size(obj.meshCoords,1), 1);
            if ~isempty(obj.dirichletDofs)
                dirNodes = ceil(obj.dirichletDofs / obj.dofsPerNode);
                isDirNode(dirNodes) = true;
            end
            
            for i = 1:nSub
                coords = obj.localMeshes{i}.coord;
                [~, globalNodes] = ismembertol(coords, obj.meshCoords, tol, 'ByRows', true);
                
                minX = min(coords(:,1)); maxX = max(coords(:,1));
                minY = min(coords(:,2)); maxY = max(coords(:,2));
                
                % Identify corner nodes (Currently the only Primal DoFs)
                isBL = abs(coords(:,1)-minX) < tol & abs(coords(:,2)-minY) < tol;
                isBR = abs(coords(:,1)-maxX) < tol & abs(coords(:,2)-minY) < tol;
                isTL = abs(coords(:,1)-minX) < tol & abs(coords(:,2)-maxY) < tol;
                isTR = abs(coords(:,1)-maxX) < tol & abs(coords(:,2)-maxY) < tol;
                
                localPrimalNodes = find(isBL | isBR | isTL | isTR);
                
                % Identify Dual (shared/Dirichlet) and Remaining nodes
                isShared       = nodeMultiplicity(globalNodes) > 1;
                isDir          = isDirNode(globalNodes);
                localDualNodes = setdiff(find(isShared | isDir), localPrimalNodes);
                localRemNodes  = setdiff(1:size(coords,1), [localPrimalNodes; localDualNodes]);
                
                % Map to DoFs
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
            % Converts physical node indices to global DoF indices
            nodes = nodes(:);
            dofs  = zeros(length(nodes) * obj.dofsPerNode, 1);
            for d = 1:obj.dofsPerNode
                dofs(d:obj.dofsPerNode:end) = (nodes - 1) * obj.dofsPerNode + d;
            end
        end
        
        function [Krr, Krp, Kpr, Kpp, fR, fP, Tdr] = splitLocalMatrices(obj, subId)
            % Splits the local stiffness matrix and force vector into r, p, and d sets
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
        
        function Ap = createAp(obj, subId, allPrimals)
            % Generates the Boolean assembly matrix for Primal DoFs
            pGlobal          = obj.primalDofsGlobal{subId};
            numPrimalsGlobal = length(allPrimals);
            numPrimalsLocal  = length(pGlobal);
            
            [~, rows] = ismember(pGlobal, allPrimals);
            cols = (1:numPrimalsLocal)';
            vals = ones(numPrimalsLocal, 1);
            
            Ap = sparse(rows, cols, vals, numPrimalsGlobal, numPrimalsLocal);
        end
        
        function activeIdx = getActivePrimalDofs(obj, allPrimals)
            % Excludes Dirichlet boundary conditions from the active primal set
            isDirichlet = ismember(allPrimals, obj.dirichletDofs);
            activeIdx   = find(~isDirichlet);
        end
        
        function Sdd = computeLocalSchur(obj, subId)
            % Computes the local Schur complement S_dd for the preconditioner
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