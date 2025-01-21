classdef Triangle_Constant < Interpolation
       
    methods (Access = public)

        function obj = Triangle_Constant(cParams)
            obj.init(cParams);
        end

    end

    methods (Access = protected)

        function computeParams(obj)
            obj.type = 'TRIANGLE';
            obj.ndime = 2;
            obj.nnode = 1;
            obj.pos_nodes = [1/3 1/3];
        end
        
        function shape = evaluateShapeFunctions(obj,xV)
            ngaus = size(xV,2);
            nelem = size(xV,3);
            N = ones(obj.nnode,ngaus, nelem);
            shape = N;
        end
        
        function deriv = evaluateShapeDerivatives(obj,xV)
            ngaus = size(xV,2);
            nelem = size(xV,3);
            dN = zeros(obj.ndime,obj.nnode,ngaus, nelem);
            deriv = dN;
        end

    end
end
