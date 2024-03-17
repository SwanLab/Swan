classdef Line_Constant < Interpolation

    methods (Access = public)

        function obj = Line_Constant(cParams)
            obj.init(cParams);
        end

    end

    methods (Access = protected)

        function computeParams(obj)
            obj.type = 'LINE';
            obj.ndime = 2;
            obj.nnode = 1;
            obj.pos_nodes = [0 0];
        end

        function shape = evaluateShapeFunctions(obj, xV)
            ngaus = length(xV);
            nelem = size(xV,3);
            N = ones(1,ngaus, nelem);
            shape = N;
        end

        function deriv = evaluateShapeDerivatives(obj, xV)
            ngaus = length(xV);
            nelem = size(xV,3);
            dN = zeros(2,ngaus, nelem);
            deriv = dN;
        end
        
    end
    
end