classdef tempFetiDPPoisson < handle
    % FETIDPPOISSON Class for the algebraic logic of the FETI-DP method.
    
    properties (Access = private)
        % Subdomain & Mesh Data
        localStiffness  
        localForces     
        localMeshes     
        meshCoords      
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
        dualIdxLocal       % Cell: Global indices within allDuals for each subdomain
        dualSignsLocal     % Cell: Boolean signs (+1/-1) for interface assembly
        localSchurBlocks   % Cell: Local Schur complement blocks (S_dd)
        dualWeights        % Array: Multiplicity weights for the Dirichlet preconditioner
    end
    
    % =========================================================
    % PUBLIC METHODS: INITIALIZATION & MAIN SOLVER
    % =========================================================
    methods (Access = public)
        
        function obj = tempFetiDPPoisson(globalMesh, subMeshes, stiffness, forces, tol, dirDofs)
            % Constructor: Initializes the FETI-DP solver and extracts DoFs
            obj.meshCoords     = globalMesh.coord;
            obj.localMeshes    = subMeshes;
            obj.localStiffness = stiffness;
            obj.localForces    = forces;
            obj.nodeTol        = tol;
            
            subdomainsVec      = [size(subMeshes, 2) size(subMeshes, 1)];
            obj.numSubdomains  = prod(subdomainsVec);
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
                Ap  = obj.createAp(subId, allPrimals);
                
                % Compute Boolean interface signs and local indices (Alternative to Bd)
                dGlobal       = obj.dualDofsGlobal{subId};
                [~, dualRows] = ismember(dGlobal, allDuals);
                isFirst       = ~visitedDuals(dGlobal);
                
                dualSigns           = ones(length(dGlobal), 1);
                dualSigns(~isFirst) = -1;
                visitedDuals(dGlobal) = true;
                
                % Store for reconstruction and preconditioner logic
                obj.dualIdxLocal{subId}   = dualRows;
                obj.dualSignsLocal{subId} = dualSigns;
                                
                % Assemble Primal Schur Complement
                sppLoc = Kpp - Kpr * Urp;
                SPP    = SPP + Ap * sppLoc * Ap';
                
                % Assemble Dual Interface system
                if ~isempty(dGlobal)
                    fDual(dualRows, dualRows) = fDual(dualRows, dualRows) + (dualSigns) .* (Urd'*Tdr') .* dualSigns';
                    dBar(dualRows)            = dBar(dualRows) + dualSigns .* (Urd' * fR);
                    BrKrrInvKrp(dualRows, :)  = BrKrrInvKrp(dualRows, :) + dualSigns .* (Tdr * Urp) * Ap';
                end
                
                rhsPrimal = rhsPrimal + Ap * (fP - Urp' * fR);
            end
                        
            activeDofs = obj.getActivePrimalDofs(allPrimals);   
            
            % Reduce system to active (non-Dirichlet) Primal DoFs
            SPPActive = SPP(activeDofs, activeDofs);                    
            BrActive  = BrKrrInvKrp(:, activeDofs);               
            rhsActive = rhsPrimal(activeDofs);                     
            
            % Final formulation of the dual interface matrix (fMat) and RHS (dBar)
            fMat = fDual + BrActive * (SPPActive \ BrActive');    
            dBar = dBar  - BrActive * (SPPActive \ rhsActive);
            
            % Setup matrices and weights for the Dirichlet preconditioner
            obj.setupPreconditioner();
        end
        
        function uFull = reconstructGlobalSolution(obj, lambdaSol, numNodes)
            % Reconstructs the full scalar field 'u' using the computed multipliers (lambda)
            allPrimals = unique(vertcat(obj.primalDofsGlobal{:}));
            numP       = length(allPrimals);
            SPP        = sparse(numP, numP);
            rhsPrimal  = zeros(numP, 1);
            
            UrpCell      = cell(obj.numSubdomains, 1);
            UrdCell      = cell(obj.numSubdomains, 1);
            KrrInvFrCell = cell(obj.numSubdomains, 1);
                                    
            % Step 1: Accumulate primal contributions and pre-calculate local inversions
            for subId = 1:obj.numSubdomains
                [Krr, Krp, Kpr, Kpp, fR, fP, Tdr] = obj.splitLocalMatrices(subId);
                Ap = obj.createAp(subId, allPrimals);
                
                Urp      = Krr \ Krp; 
                Urd      = Krr \ Tdr';
                KrrInvFr = Krr \ fR;
                
                UrpCell{subId}      = Urp;
                UrdCell{subId}      = Urd;
                KrrInvFrCell{subId} = KrrInvFr;
                
                sppLoc = Kpp - Kpr * Urp;
                SPP    = SPP + Ap * sppLoc * Ap';
                
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
            
            uP             = zeros(numP, 1);
            uP(activeDofs) = sppActive \ rhsActive;
            
            % Step 2: Reconstruct local Remainder DoFs (Internal + Duals)
            uFull = zeros(numNodes, 1); 
            
            for subId = 1:obj.numSubdomains
                pGlobal = obj.primalDofsGlobal{subId};
                dGlobal = obj.dualDofsGlobal{subId};
                rGlobal = obj.remDofsGlobal{subId};
                
                Ap        = obj.createAp(subId, allPrimals);
                dualRows  = obj.dualIdxLocal{subId};
                dualSigns = obj.dualSignsLocal{subId};
                lambdaLoc = dualSigns .* lambdaSol(dualRows);
      
                uPLoc    = Ap' * uP;
                uRemLoc  = KrrInvFrCell{subId} - UrpCell{subId} * uPLoc - UrdCell{subId} * lambdaLoc;
                
                % Scatter local results to the global solution vector
                uFull(pGlobal)      = uPLoc;
                allRemGlobal        = [rGlobal; dGlobal];
                uFull(allRemGlobal) = uRemLoc;
            end
        end
    end
    
    % =========================================================
    % PUBLIC METHODS: PRECONDITIONER
    % =========================================================
    methods (Access = public)
        
        function z = applyDirichletPrecond(obj, r)
            % Applies the Dirichlet preconditioner iteratively (Fast, Vectorized)
            numDuals = length(r);
            z = zeros(numDuals, 1);
            
            for subId = 1:obj.numSubdomains
                dualRows  = obj.dualIdxLocal{subId};
                if isempty(dualRows) 
                    continue; 
                end
                
                dualSigns = obj.dualSignsLocal{subId};
                tic;
                Sdd = obj.localSchurBlocks{subId};
                toc
                %[Kii, Kid, Kdi, Kdd]  = obj.computeLocalMatrices(subId); 
                w = obj.dualWeights(dualRows);   
                
                % Vectorized multiplication: equivalent to B_d * D * S_dd * D * B_d^T * r
                rLoc = dualSigns .* r(dualRows);
                
                z(dualRows) = z(dualRows) + dualSigns .* (w .* (Sdd * (w .* rLoc)));
                
                % % dualSigns .* (w .* ((Kdd - Kdi * (Kii \ Kid)) * (w .* rLoc)));
                % rLoc = dualSigns .* r(dualRows);
                % z1 = Kid * (w .* rLoc);
                % z2 = Kii \ z1;
                % z3 = Kdi * z2;
                % z4 = Kdd * (w .* rLoc) - z3;
                % z(dualRows) = z(dualRows) + dualSigns .* (w .* z4);
            end
        end
        
        function M = buildPrecondMatrix(obj)
            % Builds the explicit dense matrix M of the preconditioner for analysis
            numDuals = obj.getNumDuals();
            M        = zeros(numDuals); 
            
            for subId = 1:obj.numSubdomains
                dualRows  = obj.dualIdxLocal{subId};
                if isempty(dualRows)
                    continue; 
                end
                
                dualSigns = obj.dualSignsLocal{subId};
                Sdd       = obj.localSchurBlocks{subId};
                % [Kii, Kid, Kdi, Kdd] = obj.computeLocalMatrices(subId);
                % Sdd       = Kdd - Kdi * (Kii \ Kid);
                w         = obj.dualWeights(dualRows);   
                
                % Include boolean signs (B_d) and weights (w) in a diagonal matrix
                W_signed = diag(w .* dualSigns);
                
                % Additive assembly of matrix M
                M(dualRows, dualRows) = M(dualRows, dualRows) + W_signed * Sdd * W_signed;
            end
        end
               
        function n = getNumDuals(obj)
            % Retrieves the total number of unique global Dual DoFs
            allPrimals = unique(vertcat(obj.primalDofsGlobal{:}));
            allDuals   = setdiff(unique(vertcat(obj.dualDofsGlobal{:})), allPrimals);
            n          = length(allDuals);
        end
        
    end
    
    % =========================================================
    % PUBLIC METHODS: VISUALIZATION
    % =========================================================
    methods (Access = public)
        function visualizeFetiNodes(obj)
            % Visualizes the Domain Decomposition mapping (Primal, Dual, and Internal Nodes)
            allPrimalDofs = unique(vertcat(obj.primalDofsGlobal{:}));
            allDualDofs   = unique(vertcat(obj.dualDofsGlobal{:}));
            allRemDofs    = unique(vertcat(obj.remDofsGlobal{:}));
            
            primalNodes = allPrimalDofs; 
            dualNodes   = setdiff(allDualDofs, primalNodes);
            remNodes    = setdiff(allRemDofs, union(primalNodes, dualNodes));
            
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
            title('FETI-DP Node Distribution (Poisson)');
            xlabel('X'); ylabel('Y');
            grid on; hold off;
        end
    end
    
    % =========================================================
    % PRIVATE METHODS: SUBDOMAIN & MATH OPERATIONS
    % =========================================================
    methods (Access = private)
        
        % function extractFetiDofs(obj)
        %     % Evaluates geometries and extracts Primal, Dual, and Internal DoFs
        %     nSub    = prod(obj.numSubdomains);
        %     tol     = obj.nodeTol;
        % 
        %     pGlobal = cell(nSub,1); 
        %     dGlobal = cell(nSub,1); 
        %     rGlobal = cell(nSub,1);
        %     pLocal  = cell(nSub,1); 
        %     dLocal  = cell(nSub,1); 
        %     rLocal  = cell(nSub,1);
        % 
        %     for i = 1:nSub
        %         coords = obj.localMeshes{i}.coord;
        %         connec = obj.localMeshes{i}.connec;
        %         nNodes = size(coords, 1);
        % 
        %         % Local-to-global mapping
        %         [~, globalNodes] = ismembertol(coords, obj.meshCoords, tol, 'ByRows', true);
        % 
        %         % Extract subdomain boundary edges
        %         allEdges       = [connec(:,[1,2]); connec(:,[2,3]); connec(:,[3,1])];
        %         allEdgesSorted = sort(allEdges, 2);
        %         [~, ia, ic]    = unique(allEdgesSorted, 'rows');
        %         edgeCount      = accumarray(ic, 1);
        %         boundaryEdges  = allEdges(ia(edgeCount == 1), :);
        %         boundaryLocalIdx = unique(boundaryEdges(:));
        % 
        %        % Obtain local boundary meshes
        %         bMeshes = obj.localMeshes{i}.createBoundaryMesh();
        % 
        %         % PRIMAL: Endpoints (corners) of each boundary mesh
        %         primalGlobalSet = [];
        %         for b = 1:length(bMeshes)
        %             gc          = bMeshes{b}.globalConnec;  
        %             allN        = gc(:);
        %             uN          = unique(allN);
        %             deg         = accumarray(allN, 1, [max(allN), 1]);
        %             extremLocal = uN(deg(uN) == 1);
        %             primalGlobalSet = [primalGlobalSet; globalNodes(extremLocal)];
        %         end
        %         primalGlobalSet = unique(primalGlobalSet);
        % 
        %         % Convert to local indices
        %         [~, primalLocalIdx] = ismember(primalGlobalSet, globalNodes);
        %         primalLocalIdx      = primalLocalIdx(primalLocalIdx > 0);
        % 
        %        % DUAL: Interface/boundary nodes excluding primals
        %         bCoords = coords(boundaryLocalIdx, :);
        %         minX = min(bCoords(:,1));  maxX = max(bCoords(:,1));
        %         minY = min(bCoords(:,2));  maxY = max(bCoords(:,2));
        % 
        %         isOnExtreme = (abs(coords(:,1) - minX) < tol) | ...
        %                       (abs(coords(:,1) - maxX) < tol) | ...
        %                       (abs(coords(:,2) - minY) < tol) | ...
        %                       (abs(coords(:,2) - maxY) < tol);
        % 
        %         extremeLocalIdx = find(isOnExtreme);
        %         dualLocalIdx    = setdiff(intersect(extremeLocalIdx, boundaryLocalIdx), primalLocalIdx);
        % 
        %         % REMAINING: Internal nodes
        %         remLocalIdx = setdiff((1:nNodes)', [primalLocalIdx; dualLocalIdx]);
        % 
        %         pGlobal{i} = globalNodes(primalLocalIdx);
        %         dGlobal{i} = globalNodes(dualLocalIdx);
        %         rGlobal{i} = globalNodes(remLocalIdx);
        % 
        %         pLocal{i}  = primalLocalIdx;
        %         dLocal{i}  = dualLocalIdx;
        %         rLocal{i}  = remLocalIdx;
        %     end
        % 
        %     obj.primalDofsGlobal = pGlobal; obj.primalDofsLocal = pLocal;
        %     obj.dualDofsGlobal   = dGlobal; obj.dualDofsLocal   = dLocal;
        %     obj.remDofsGlobal    = rGlobal; obj.remDofsLocal    = rLocal;
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
        
                % Local -> global map
                [~,globalNodes] = ismembertol( ...
                    coords, obj.meshCoords, tol, 'ByRows', true);
        
                % Boundary meshes
                bMeshes = obj.localMeshes{i}.createBoundaryMesh();
        
                % All boundary nodes
                boundaryNodes = [];
        
                % Corner nodes (primal)
                primalNodes = [];
        
                for b = 1:length(bMeshes)
        
                    gc = bMeshes{b}.globalConnec;
        
                    boundaryNodes = [boundaryNodes; gc(:)];
        
                    % Endpoints of 1D boundary chain
                    uN  = unique(gc(:));
                    deg = accumarray(gc(:),1);
        
                    corners = uN(deg(uN)==1);
        
                    primalNodes = [primalNodes; corners];
        
                end
        
                boundaryNodes = unique(boundaryNodes);
                primalNodes   = unique(primalNodes);
        
                % Dual = boundary minus primal
                dualNodes = setdiff(boundaryNodes, primalNodes);
        
                % Remaining = internal
                remNodes = setdiff((1:nNodes)', ...
                            [primalNodes; dualNodes]);
        
                % Store local
                pLocal{i} = primalNodes;
                dLocal{i} = dualNodes;
                rLocal{i} = remNodes;
        
                % Store global
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

        function [Krr, Krp, Kpr, Kpp, fR, fP, Tdr] = splitLocalMatrices(obj, subId)
            % Splits the local stiffness matrix and force vector into r, p, and d sets
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
            
            % Create Trace operator boolean matrix (T_dr)
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
            % Generates the Boolean assembly matrix for Primal DoFs
            pGlobal = obj.primalDofsGlobal{subId};
            numPrimalsGlobal = length(allPrimals);
            numPrimalsLocal  = length(pGlobal);
            
            [~, rows] = ismember(pGlobal, allPrimals);
            cols      = (1:numPrimalsLocal)';
            vals      = ones(numPrimalsLocal, 1);
            
            ap = sparse(rows, cols, vals, numPrimalsGlobal, numPrimalsLocal);
        end
        
        function activeIdx = getActivePrimalDofs(obj, allPrimals)
            % Excludes Dirichlet boundary conditions from the active primal set
            isDirichlet = ismember(allPrimals, obj.dirichletDofs);
            activeIdx   = find(~isDirichlet);  
        end
        
        function Sdd = computeLocalSchur(obj, subId)
            % Computes the local Schur complement S_dd for the preconditioner
            kMat = obj.localStiffness{subId};
            rLoc = obj.remDofsLocal{subId};   % Internal
            dLoc = obj.dualDofsLocal{subId};  % Dual Interface

            Kii = kMat(rLoc, rLoc);
            Kid = kMat(rLoc, dLoc);
            Kdi = kMat(dLoc, rLoc);
            Kdd = kMat(dLoc, dLoc);

            Sdd = Kdd - Kdi * (Kii \ Kid);
        end

        function [Kii, Kid, Kdi, Kdd] = computeLocalMatrices(obj, subId)
            % Computes the local Schur complement S_dd for the preconditioner
            kMat = obj.localStiffness{subId};
            rLoc = obj.remDofsLocal{subId};   % Internal
            dLoc = obj.dualDofsLocal{subId};  % Dual Interface
        
            Kii = kMat(rLoc, rLoc);
            Kid = kMat(rLoc, dLoc);
            Kdi = kMat(dLoc, rLoc);
            Kdd = kMat(dLoc, dLoc);
        
            %Sdd = Kdd - Kdi * (Kii \ Kid);
        end
        
        function setupPreconditioner(obj)
            % Computes preconditioner weights and pre-assembles Schur blocks
            numDuals     = obj.getNumDuals();
            multiplicity = zeros(numDuals, 1);
            
            obj.localSchurBlocks = cell(obj.numSubdomains, 1);
            
            for subId = 1:obj.numSubdomains
                obj.localSchurBlocks{subId} = obj.computeLocalSchur(subId);
                
                dualRows               = obj.dualIdxLocal{subId};
                multiplicity(dualRows) = multiplicity(dualRows) + 1;
            end
            
            obj.dualWeights = 1 ./ multiplicity;
        end
        
    end
end