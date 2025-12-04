classdef LineDegree8 < Interpolation
    % 1D 8th-degree Lagrange interpolation using Chebyshev nodes in [-1,1]
    
    methods (Access = public)
        function obj = LineDegree8(cParams)
            obj.init(cParams);
            obj.computeParams();
        end
    end
    
    methods (Access = protected)
        
        function computeParams(obj)
            obj.ndime = 1;
            obj.nnode = 9; % degree 8 => 9 nodes
            n = obj.nnode;

            % Chebyshev nodes in [-1,1]
            k = (1:n)';
            obj.pos_nodes = -cos((2*k - 1)*pi/(2*n));
        end

        function shape = evaluateShapeFunctions(obj,posgp)
            % Vectorized, MATLAB-style (like your quadratic version)
            ngaus = size(posgp,2);
            nelem = size(posgp,3);
            s = posgp(1,:,:); % size [1,ngauss,nelem]
            shape = zeros(obj.nnode,ngaus,nelem);

            xnodes = obj.pos_nodes(:);
            n = obj.nnode;

            % Precompute denominators
            denom = ones(n,1);
            for i = 1:n
                denom(i) = prod(xnodes(i) - xnodes([1:i-1,i+1:end]));
            end

            % Compute each Lagrange polynomial in vectorized form
            for i = 1:n
                others = xnodes([1:i-1,i+1:end]);
                % Elementwise product over (s - others(j))
                term = ones(size(s));
                for j = 1:numel(others)
                    term = term .* (s - others(j));
                end
                shape(i,:,:) = term / denom(i);
            end
        end

        function deriv = evaluateShapeDerivatives(obj,posgp)
            ngaus = size(posgp,2);
            nelem = size(posgp,3);
            s = posgp(1,:,:);
            n = obj.nnode;
            xnodes = obj.pos_nodes(:);

            deriv = zeros(obj.ndime,n,ngaus,nelem);

            % Precompute denominators
            denom = ones(n,1);
            for i = 1:n
                denom(i) = prod(xnodes(i) - xnodes([1:i-1,i+1:end]));
            end

            % Derivative of Lagrange polynomials
            for i = 1:n
                dLi = zeros(size(s));
                for k = 1:n
                    if k == i, continue; end
                    others = xnodes([1:i-1,i+1:end]);
                    others(others == xnodes(k)) = []; % exclude k-th node
                    term = ones(size(s));
                    for j = 1:numel(others)
                        term = term .* (s - others(j));
                    end
                    dLi = dLi + term / (xnodes(i) - xnodes(k));
                end
                deriv(1,i,:,:) = dLi / denom(i);
            end
        end
    end
end
