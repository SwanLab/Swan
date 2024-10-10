classdef Triangle_RaviartThomas < Interpolation

    methods (Access = public)

        function obj = Triangle_RaviartThomas(cParams)
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
            shape = zeros(obj.nnode,ngaus,nelem,2);
            s = xV(1,:,:);
            t = xV(2,:,:);
            I = ones(size(s));

            shape(1,:,:,1) = s;
            shape(1,:,:,2) = t-I;
            shape(2,:,:,1) = s;
            shape(2,:,:,2) = t;
            shape(3,:,:,1) = s-I;
            shape(3,:,:,2) = t;
        end

        function deriv = evaluateShapeDerivatives(obj, xV)
            ngaus = size(xV,2);
            nelem = size(xV,3);
            deriv = zeros(obj.ndime,obj.nnode,ngaus,nelem,2);
            s = xV(1,:,:);
            I = ones(size(s));
            O = zeros(size(s));

            deriv(1,1,:,:,1) = I;
            deriv(1,1,:,:,2) = O;
            deriv(1,2,:,:,1) = I;
            deriv(1,2,:,:,2) = O;
            deriv(1,3,:,:,1) = I;
            deriv(1,3,:,:,2) = O;

            deriv(2,1,:,:,1) = O;
            deriv(2,1,:,:,2) = I;
            deriv(2,2,:,:,1) = O;
            deriv(2,2,:,:,2) = I;
            deriv(2,3,:,:,1) = O;
            deriv(2,3,:,:,2) = I;
        end

    end

end
