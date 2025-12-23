classdef LineDegree8 < Interpolation
    
    properties (Access = protected)
        % Matrix of coefficients for the Derivatives [9 x 8]
        % Row i contains the polynomial coefficients for the derivative of shape function i
        deriv_coeffs 
        
        % Matrix of coefficients for the Shape Functions [9 x 9] (Optional, for speed)
        shape_coeffs
    end
    
    methods (Access = public)
        function obj = LineDegree8(cParams)
            obj.init(cParams);
            obj.computeParams();
        end
    end
    
    methods (Access = protected)
        
        function computeParams(obj)
            obj.ndime = 1;
            obj.nnode = 9;
            n = obj.nnode;
            degree = n - 1; % Degree 8
            
            % 1. Chebyshev-Gauss-Lobatto nodes
            k = (1:n)';
            x = -cos(pi * (k-1) / degree);
            obj.pos_nodes = x;
            
            % 2. Pre-compute Polynomial Coefficients (The "Monomial Basis")
            % We find the coefficients that define each Lagrange polynomial.
            obj.shape_coeffs = zeros(n, n);      % Degree 8 (9 coeffs)
            obj.deriv_coeffs = zeros(n, degree); % Derivative is Degree 7 (8 coeffs)
            
            for i = 1:n
                % Create the "target": 1 at node i, 0 at all other nodes
                y = zeros(size(x));
                y(i) = 1;
                
                % Fit polynomial (find coefficients for this shape function)
                p = polyfit(x, y, degree); 
                
                % Store Shape Coeffs
                obj.shape_coeffs(i, :) = p;
                
                % Calculate and Store Derivative Coeffs
                dp = polyder(p);
                obj.deriv_coeffs(i, :) = dp;
            end
        end

        function shape = evaluateShapeFunctions(obj, posgp)
            % s: [1, nGaus, nElem]
            s = posgp(1,:,:);
            
            % Flatten s for matrix operation
            s_flat = s(:);
            M = length(s_flat);
            
            % 1. Create Powers Matrix (Vandermonde-like)
            % We need [s^8, s^7, ... s^1, 1]
            % Doing this explicitly is faster than `vander`
            pow_matrix = zeros(9, M); 
            pow_matrix(9, :) = 1;      % s^0
            pow_matrix(8, :) = s_flat; % s^1
            
            % Compute higher powers iteratively (faster than .^)
            for p = 2:8
                pow_matrix(9-p, :) = s_flat .* pow_matrix(9-p+1, :).'; 
            end
            
            % 2. Matrix Multiply: [9x9] * [9xM] = [9xM]
            shapes_flat = obj.shape_coeffs * pow_matrix;
            
            % 3. Reshape result
            shape = reshape(shapes_flat, [obj.nnode, size(s,2), size(s,3)]);
        end

        function deriv = evaluateShapeDerivatives(obj, posgp)
            % s: [1, nGaus, nElem]
            s = posgp(1,:,:);
            
            % Flatten s
            s_flat = s(:);
            M = length(s_flat);
            
            % 1. Create Powers Matrix for Derivative (Degree 7)
            % We need [s^7, s^6, ..., s^1, 1]
            pow_matrix = zeros(8, M);
            pow_matrix(8, :) = 1;      % s^0
            pow_matrix(7, :) = s_flat; % s^1
            
            % Compute higher powers iteratively
            for p = 2:7
                pow_matrix(8-p, :) = s_flat .* pow_matrix(8-p+1, :).';
            end
            
            % 2. Matrix Multiply: [9x8] * [8xM] = [9xM]
            % This computes all derivatives for all nodes at all points instantly
            d_shapes_flat = obj.deriv_coeffs * pow_matrix;
            
            % 3. Reshape result [nDime, nNode, nGauss, nElem]
            deriv = reshape(d_shapes_flat, [1, obj.nnode, size(s,2), size(s,3)]);
        end
    end
end