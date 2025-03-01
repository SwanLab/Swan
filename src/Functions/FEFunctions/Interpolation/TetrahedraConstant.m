classdef TetrahedraConstant < Interpolation

    methods (Access = public)

        function obj = TetrahedraConstant(cParams)
            obj.init(cParams);
        end

    end

    methods (Access = protected)

        function computeParams(obj)
            obj.type = 'TETRAHEDRA';
            obj.ndime = 3;
            obj.nnode = 1;
            obj.pos_nodes = [1/4 1/4 1/4];
        end

        function shape = evaluateShapeFunctions(obj,xV)
            ngaus = size(xV,2);
            nelem = size(xV,3);
            N = ones(obj.nnode,ngaus,nelem);
            shape = N;
        end

        function deriv = evaluateShapeDerivatives(obj,xV)
            ngaus = size(xV,2);
            nelem = size(xV,3);
            dN = zeros(obj.ndime,obj.nnode,ngaus,nelem);
            deriv = dN;
        end

    end
end
