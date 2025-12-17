classdef LineDegree8 < Interpolation
    
    % Add a property to store denominators so we compute them only once
    properties (Access = protected)
        denom
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
            
            % 1. Chebyshev-Gauss-Lobatto nodes
            k = (1:n)';
            obj.pos_nodes = -cos(pi * (k-1) / (n-1));
            
            % 2. Precompute Denominators (Optimization)
            % This is O(N^2) but runs only once during initialization
            x = obj.pos_nodes;
            obj.denom = zeros(n,1);
            for i = 1:n
                % Product of (xi - xj) for all j != i
                obj.denom(i) = prod(x(i) - x([1:i-1, i+1:end]));
            end
        end

        function shape = evaluateShapeFunctions(obj, posgp)
            % s: [1, nGaus, nElem]
            s = posgp(1,:,:); 
            x = obj.pos_nodes; % [nNode, 1]
            n = obj.nnode;
            
            % 1. Create Difference Matrix [nNode, nGaus, nElem]
            % Uses Implicit Expansion (s - x) automatically broadcasts dimensions
            Diff = s - x; 
            
            shape = zeros(n, size(s,2), size(s,3));

            % 2. Compute Shape Functions
            for i = 1:n
                % Vectorized Product: Compute product across the node dimension (dim 1)
                % We select all rows of Diff except row 'i'
                numerator = prod(Diff([1:i-1, i+1:end], :, :), 1);
                
                % Apply precomputed denominator
                shape(i,:,:) = numerator / obj.denom(i);
            end
        end

        function deriv = evaluateShapeDerivatives(obj, posgp)
            s = posgp(1,:,:);
            x = obj.pos_nodes;
            n = obj.nnode;
            
            % 1. Difference Matrix [nNode, nGaus, nElem]
            Diff = s - x;
            
            deriv = zeros(obj.ndime, n, size(s,2), size(s,3));

            % 2. Compute Derivatives (Sum of Products)
            for i = 1:n
                % We need the sum of terms where each term excludes 'i' and one 'k'
                sum_terms = zeros(size(s));
                
                % Create indices for row 'i' to exclude it easily below
                all_idx = 1:n;
                idx_no_i = all_idx(all_idx ~= i); 
                
                for k_idx = 1:length(idx_no_i)
                    k = idx_no_i(k_idx);
                    
                    % We want product of (s - xj) for all j != i AND j != k
                    % This creates the indices excluding i and k
                    idx_prod = idx_no_i(idx_no_i ~= k);
                    
                    % Vectorized product over the remaining nodes
                    term = prod(Diff(idx_prod, :, :), 1);
                    sum_terms = sum_terms + term;
                end
                
                deriv(1,i,:,:) = sum_terms / obj.denom(i);
            end
        end
    end
end