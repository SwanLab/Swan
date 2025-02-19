classdef LineConstant < Interpolation

    methods (Access = public)

        function obj = LineConstant(cParams)
            obj.init(cParams);
        end

    end

    methods (Access = protected)

        function computeParams(obj)
            obj.ndime = 1;
            obj.nnode = 1;
            obj.pos_nodes = [0 0];
        end

        function shape = evaluateShapeFunctions(obj,posgp)
            ngaus = length(posgp);
            nelem = size(posgp,3);
            N = ones(obj.nnode,ngaus, nelem);
            shape = N;
        end

        function deriv = evaluateShapeDerivatives(obj,posgp)
            ngaus = length(posgp);
            nelem = size(posgp,3);
            dN = zeros(obj.ndime,ngaus, nelem);
            deriv = dN;
        end

    end

end