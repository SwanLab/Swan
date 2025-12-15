classdef LineDegree8V2 < Interpolation
    
    methods (Access = public)
        function obj = LineDegree8V2(cParams)
            obj.init(cParams);
            obj.computeParams();
        end
    end
    
    methods (Access = protected)
        
        function computeParams(obj)
            obj.ndime = 1;
            obj.nnode = 9;
            n = obj.nnode;
            k = (1:n)';
            obj.pos_nodes = -cos(pi * (k-1) / (n-1));
        end

        function shape = evaluateShapeFunctions(obj,posgp)
            ngaus = size(posgp,2);
            nelem = size(posgp,3);
            s = posgp(1,:,:); 
            shape = zeros(obj.nnode,ngaus,nelem);

            xnodes = obj.pos_nodes(:);
            n = obj.nnode;

            denom = ones(n,1);
            for i = 1:n
                denom(i) = prod(xnodes(i) - xnodes([1:i-1,i+1:end]));
            end

            for i = 1:n
                others = xnodes([1:i-1,i+1:end]);
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

            denom = ones(n,1);
            for i = 1:n
                denom(i) = prod(xnodes(i) - xnodes([1:i-1,i+1:end]));
            end

            for i = 1:n
                dLi = zeros(size(s));
                for k = 1:n
                    if k == i, continue; end
                    
                    others = xnodes([1:i-1,i+1:end]);
                    others(others == xnodes(k)) = [];
                    
                    term = ones(size(s));
                    for j = 1:numel(others)
                        term = term .* (s - others(j));
                    end
                    dLi = dLi + term;
                end
                deriv(1,i,:,:) = dLi / denom(i);
            end
        end
    end
end