classdef Line_Constant < Interpolation

    methods (Access = public)

        function obj = Line_Constant(cParams)
            obj.init(cParams);
            obj.computeParams();
        end
        
        function shape = computeShapeFunctions(obj,posgp)
            ngaus = length(posgp);
            N = ones(1,ngaus);
            shape = N;
        end
        
        function deriv = computeShapeDerivatives(obj,posgp)
            ngaus = length(posgp);
            dN = zeros(2,ngaus);
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