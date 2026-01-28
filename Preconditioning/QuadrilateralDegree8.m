classdef QuadrilateralDegree8 < Interpolation

    properties (Access = protected)
        % We keep the 1D coefficients because the 2D basis is a tensor product
        % of the 1D basis.
        % deriv_coeffs: [9 x 8]
        % shape_coeffs: [9 x 9]
        deriv_coeffs
        shape_coeffs
        perm_map
    end

    methods (Access = public)
        function obj = QuadrilateralDegree8(cParams)
            obj.init(cParams);
            obj.get_permutation_p8();
            obj.computeParams();
        end
    end

    methods (Access = protected)

        function computeParams(obj)
            obj.ndime = 2;
            obj.nnode = 81;
            n1D = 9;
            degree = n1D - 1;

            % 1. Chebyshev-Gauss-Lobatto nodes (1D)
            k = (1:n1D)';
            x1D = -cos(pi * (k-1) / degree);

            % 2. Create 2D Tensor Product Grid
            [X, Y] = meshgrid(x1D, x1D);

            % FIX: Transpose X and Y so that when we flatten with (:),
            % we traverse the ROW (x) first, then jump to the next column (y).
            X = X';
            Y = Y';

            % Now X(:) changes fastest (x1, x2, x3...), matching your loop logic
            obj.pos_nodes = [X(:), Y(:)];

            % 3. Pre-compute 1D Polynomial Coefficients
            % (Exactly the same logic as LineDegree8)
            obj.shape_coeffs = zeros(n1D, n1D);
            obj.deriv_coeffs = zeros(n1D, degree);

            for i = 1:n1D
                % Target: 1 at node i, 0 at others
                y = zeros(size(x1D));
                y(i) = 1;

                % Fit polynomial
                p = polyfit(x1D, y, degree);
                obj.shape_coeffs(i, :) = p;

                % Derivative coeffs
                dp = polyder(p);
                obj.deriv_coeffs(i, :) = dp;
            end
        end

        function shape = evaluateShapeFunctions(obj, posgp)
            % posgp: [2, nGauss, nElem]
            s = posgp(1,:,:);
            t = posgp(2,:,:);

            % 1. Compute 1D basis for s and t separately
            % We only need shapes (flag=0)
            [Ns, ~] = obj.compute1DBasis(s(:).'); % [9 x M]
            [Nt, ~] = obj.compute1DBasis(t(:).'); % [9 x M]

            M = size(Ns, 2);

            % 1. Compute Tensor Product (Generates Lexicographical Order)
            % (Use the Corrected S-Fastest version I gave you previously)
            shapes_lex = zeros(obj.nnode, M);
            for j = 1:9
                startIdx = (j-1)*9 + 1;
                endIdx   = j*9;
                shapes_lex(startIdx:endIdx, :) = Ns .* Nt(j, :);
            end

            % 2. PERMUTE to match Mesh Order (Corners -> Edges -> Internal)
            % We assume 'obj.perm_map' is precomputed using the logic above
            % shape_lex(k, :) moves to row p_map(k) ??
            % Actually, we usually want: Output(k) = Input( Map(k) )

            % If p_map lists the Grid Index for the k-th Mesh Node:
            shape = shapes_lex(obj.perm_map, :);

            % Reshape if needed
            shape = reshape(shape, [obj.nnode, size(s,2), size(s,3)]);
        end

        function deriv = evaluateShapeDerivatives(obj, posgp)
            % posgp: [2, nGauss, nElem]
            s = posgp(1,:,:);
            t = posgp(2,:,:);

            % 1. Compute 1D basis AND derivatives for s and t
            [Ns, dNs] = obj.compute1DBasis(s(:).'); % [9 x M]
            [Nt, dNt] = obj.compute1DBasis(t(:).'); % [9 x M]

            M = size(Ns, 2);
            nNodes = obj.nnode;

            % 2. Initialize Output [nDime, nNode, M]
            % We flatten nGauss x nElem into M for calculation
            % 1. Compute Tensor Product Derivatives (Lexicographical X-Fastest)
            d_shapes_lex = zeros(2, nNodes, M);

            for j = 1:9
                startIdx = (j-1)*9 + 1;
                endIdx   = j*9;

                % d/dxi component (s) -> Row 1
                d_shapes_lex(1, startIdx:endIdx, :) = dNs .* Nt(j, :);

                % d/deta component (t) -> Row 2
                d_shapes_lex(2, startIdx:endIdx, :) = Ns .* dNt(j, :);
            end

            % 2. CRITICAL FIX: Apply Permutation to Derivatives
            % We must permute the 'nNode' dimension (dimension 2)
            deriv_permuted = d_shapes_lex(:, obj.perm_map, :);

            % 3. Reshape
            deriv = reshape(deriv_permuted, [2, nNodes, size(s,2), size(s,3)]);
        end
    end

    methods (Access = private)

        function [N, dN] = compute1DBasis(obj, coord_flat)
            % Helper to compute 1D shapes and derivatives for a list of coordinates
            % coord_flat: [1 x M] Row Vector

            M = length(coord_flat);

            % --- Shape Functions (Degree 8) ---
            % Powers: [s^8, s^7, ... 1]
            pow_mat = zeros(9, M);
            pow_mat(9, :) = 1;
            pow_mat(8, :) = coord_flat;

            for p = 2:8
                pow_mat(9-p, :) = coord_flat .* pow_mat(9-p+1, :);
            end

            N = obj.shape_coeffs * pow_mat; % [9 x 9] * [9 x M] = [9 x M]

            % --- Derivatives (Degree 7) ---
            % Powers: [s^7, s^6, ... 1]
            pow_mat_d = zeros(8, M);
            pow_mat_d(8, :) = 1;
            pow_mat_d(7, :) = coord_flat;

            for p = 2:7
                pow_mat_d(8-p, :) = coord_flat .* pow_mat_d(8-p+1, :);
            end

            dN = obj.deriv_coeffs * pow_mat_d; % [9 x 8] * [8 x M] = [9 x M]
        end

        function get_permutation_p8(obj)
            % p = 8 implies 9 nodes per side. Total 81 nodes.
            % We want to map FROM Grid (i,j) TO Mesh Index.

            N = 9; % Nodes per side
            grid_indices = reshape(1:81, N, N); % 9x9 grid of indices (Lexicographical)

            % --- 1. Corners (Counter-Clockwise usually) ---
            % standard: Bottom-Left, Bottom-Right, Top-Right, Top-Left
            corners = [grid_indices(1,1), grid_indices(N,1), ...
                grid_indices(N,N), grid_indices(1,N)];

            % --- 2. Edges (Counter-Clockwise) ---
            % Bottom Edge (excluding corners)
            edge_bottom = grid_indices(2:N-1, 1);
            % Right Edge
            edge_right  = grid_indices(N, 2:N-1);
            % Top Edge (often reversed to maintain CCW loop)
            edge_top    = grid_indices(N-1:-1:2, N);
            % Left Edge
            edge_left   = grid_indices(1, N-1:-1:2);

            edges = [edge_bottom; edge_right(:); edge_top(:); edge_left(:)];

            % --- 3. Interior ---
            % The inner 7x7 block
            interior = grid_indices(2:N-1, 2:N-1);
            interior = interior(:); % Vectorize

            % Combine them
            % This vector lists the Grid Indices in the order the Mesh wants them
            obj.perm_map = [corners(:); edges(:); interior(:)];
        end
    end
end