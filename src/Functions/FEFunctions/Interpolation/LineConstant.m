classdef LineConstant < Interpolation

    methods (Access = public)

        function obj = LineConstant(cParams)
            obj.init(cParams);
            obj.computeParams();
        end
        
        function shape = computeShapeFunctions(obj,posgp)
            ngaus = length(posgp);
            nelem = size(posgp,3);
            N = ones(obj.nnode,ngaus, nelem);
            shape = N;
        end
        
        function deriv = computeShapeDerivatives(obj,posgp)
            ngaus = length(posgp);
            nelem = size(posgp,3);
            dN = zeros(obj.ndime,ngaus, nelem);
            deriv = dN;
        end

    end

    methods (Access = private)

        function computeParams(obj)
            obj.type = 'LINE';
            obj.ndime = 2;
            obj.nnode = 1;
            obj.pos_nodes = [0 0];
        end
        
    end
    
end