classdef LineDegree8 < Interpolation
    
    properties (Access = protected)
        % Matrix of derivative coefficients [nNode x (Degree)]
        % Row i contains coefficients for the derivative of shape function i
        deriv_coeffs 
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
            degree = n - 1;
            
            % 1. Chebyshev-Gauss-Lobatto nodes
            k = (1:n)';
            x = -cos(pi * (k-1) / degree);
            obj.pos_nodes = x;
            
            % 2. Pre-compute Derivative Coefficients
            % We fit a polynomial for each node (delta function) and differentiate it
            obj.deriv_coeffs = zeros(n, degree); % Derivative is degree 7 (N-2)
            
            for i = 1:n
                % Create target values (1 at node i, 0 elsewhere)
                y = zeros(size(x));
                y(i) = 1;
                
                % Get polynomial coeffs (p) for shape function i
                % and then the derivative coeffs (dp)
                p = polyfit(x, y, degree); 
                dp = polyder(p);
                
                % Store in matrix (pad if necessary, though degree is fixed)
                obj.deriv_coeffs(i, :) = dp;
            end
        end

        function shape = evaluateShapeFunctions(obj, posgp)
            % This can also be optimized similarly, but let's stick 
            % to the derivative request first.
            % (Keeping your existing logic or the previous optimization here)
            
            % ... [Previous implementation of evaluateShapeFunctions] ... 
            % For consistency, I will paste the fast Vectorized version here:
            s = posgp(1,:,:); 
            x = obj.pos_nodes; 
            n = obj.nnode;
            Diff = s - x; 
            
            % Re-calculate denom since we aren't using polyval here for consistency
            denom = zeros(n,1);
            for i = 1:n, denom(i) = prod(x(i) - x([1:i-1, i+1:end])); end

            shape = zeros(n, size(s,2), size(s,3));
            for i = 1:n
                numerator = prod(Diff([1:i-1, i+1:end], :, :), 1);
                shape(i,:,:) = numerator / denom(i);
            end
        end

        function deriv = evaluateShapeDerivatives(obj, posgp)
            % s: [1, nGaus, nElem]
            s = posgp(1,:,:); 
            
            % Flatten s for vectorized evaluation: [M x 1]
            s_flat = s(:);
            M = length(s_flat);
            
            % 1. Construct Vandermonde-like matrix for the derivative powers
            % The derivative of degree 8 is degree 7.
            % Powers needed: [s^7, s^6, ..., s^1, 1]
            current_degree = obj.nnode - 2; % 9 nodes -> degree 8 -> deriv degree 7
            
            % Create powers matrix: [M x 8]
            % This is much faster than loops
            S_powers = ones(M, current_degree + 1);
            for p = 1:current_degree
                S_powers(:, end - p) = s_flat.^p;
            end
            
            % 2. Matrix Multiply: [nNode x 8] * [8 x M] = [nNode x M]
            % deriv_coeffs is [9 x 8]
            % S_powers'    is [8 x M]
            d_shapes_flat = obj.deriv_coeffs * S_powers';
            
            % 3. Reshape back to [1, nNode, nGaus, nElem]
            deriv = reshape(d_shapes_flat, [1, obj.nnode, size(s,2), size(s,3)]);
        end
    end
end