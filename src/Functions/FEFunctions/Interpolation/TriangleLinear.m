classdef TriangleLinear < Interpolation

    methods (Access = public)

        function obj = TriangleLinear(cParams)
            obj.init(cParams);
        end

    end

    methods (Access = protected)

        function computeParams(obj)
            obj.ndime     = 2;
            obj.nnode     = 3;
            obj.pos_nodes = [0 0; 1 0; 0 1];
        end

        function shape = evaluateShapeFunctions(obj, xV)
            ngaus = size(xV,2);
            nelem = size(xV,3);
            s = xV(1,:,:);
            t = xV(2,:,:);
            I = ones(size(t));
            shape = zeros(obj.nnode,ngaus,nelem);
            shape(1,:,:) = I-s-t;
            shape(2,:,:) = s;
            shape(3,:,:) = t;
        end

        function deriv = evaluateShapeDerivatives(obj, xV)
            ngaus = size(xV,2);
            nelem = size(xV,3);
            deriv = zeros(obj.ndime,obj.nnode,ngaus,nelem);
            deriv(1,1,:,:) = -1;
            deriv(1,2,:,:) = 1;
            deriv(1,3,:,:) = 0;
            deriv(2,1,:,:) = -1;
            deriv(2,2,:,:) = 0;
            deriv(2,3,:,:) = 1;
        end

    end

end
