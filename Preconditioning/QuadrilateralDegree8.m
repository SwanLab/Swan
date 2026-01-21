classdef QuadrilateralDegree8 < Interpolation
    
    properties (Access = protected)
        % We keep the 1D coefficients because the 2D basis is a tensor product
        % of the 1D basis.
        % deriv_coeffs: [9 x 8]
        % shape_coeffs: [9 x 9]
        deriv_coeffs 
        shape_coeffs
    end
    
    methods (Access = public)
        function obj = QuadrilateralDegree8(cParams)
            obj.init(cParams);
            obj.computeParams();
        end
    end
    
    methods (Access = protected)
        
        function computeParams(obj)
            obj.ndime = 2;
            obj.nnode = 81; % 9 * 9 nodes
            n1D = 9;
            degree = n1D - 1;
            
            % 1. Chebyshev-Gauss-Lobatto nodes (1D)
            k = (1:n1D)';
            x1D = -cos(pi * (k-1) / degree);
            
            % 2. Create 2D Tensor Product Grid
            % We order them lexicographically: Loop y, then x
            [X, Y] = meshgrid(x1D, x1D);
            % Transpose to match standard FEM node ordering (x varies fastest)
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
            
            % 2. Compute Tensor Product: N = Ns(i) * Nt(j)
            % Result size: [81, M]
            % We utilize implicit expansion or looping for clarity/speed
            shapes_flat = zeros(obj.nnode, M);
            
            % Ordering: x (s) varies fastest (inner loop), y (t) varies slowest
            for j = 1:9
                % Broadcast Nt(j) across all Ns rows
                % This creates the block of shape functions for row j of the grid
                startIdx = (j-1)*9 + 1;
                endIdx   = j*9;
                shapes_flat(startIdx:endIdx, :) = Ns .* Nt(j, :);
            end
            
            % 3. Reshape
            shape = reshape(shapes_flat, [obj.nnode, size(s,2), size(s,3)]);
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
            d_shapes_flat = zeros(2, nNodes, M);
            
            % 3. Apply Tensor Product Chain Rule
            % dN/dxi  = dNs(i) * Nt(j)
            % dN/deta = Ns(i)  * dNt(j)
            
            for j = 1:9
                startIdx = (j-1)*9 + 1;
                endIdx   = j*9;
                
                % d/dxi component (s)
                d_shapes_flat(1, startIdx:endIdx, :) = dNs .* Nt(j, :);
                
                % d/deta component (t)
                d_shapes_flat(2, startIdx:endIdx, :) = Ns .* dNt(j, :);
            end
            
            % 4. Reshape
            deriv = reshape(d_shapes_flat, [2, nNodes, size(s,2), size(s,3)]);
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
    end
end