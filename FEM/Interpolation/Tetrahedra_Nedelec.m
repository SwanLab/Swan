classdef Tetrahedra_Nedelec < Interpolation

    methods (Access = public)

        function obj = Tetrahedra_Nedelec(cParams)
            obj.init(cParams);
        end

    end

    methods (Access = protected)

        function computeParams(obj)
            obj.ndime     = 3;
            obj.nnode     = 6;
            obj.pos_nodes = [0 0 0; 1 0 0; 0 1 0; 0 0 1];
        end

        function shape = evaluateShapeFunctions(obj, xV)
            ngaus = size(xV,2);
            nelem = size(xV,3);
            shape = zeros(obj.nnode,ngaus,nelem,3);
            s = xV(1,:,:);
            t = xV(2,:,:);
            u = xV(3,:,:);
            I = ones(size(s));
            O = zeros(size(s));

            shape(1,:,:,1) = I-t-u;
            shape(1,:,:,2) = s;
            shape(1,:,:,3) = s;

            shape(2,:,:,1) = t;
            shape(2,:,:,2) = I-s-u;
            shape(2,:,:,3) = t;

            shape(3,:,:,1) = u;
            shape(3,:,:,2) = u;
            shape(3,:,:,3) = I-s-t;

            shape(4,:,:,1) = t;
            shape(4,:,:,2) = -s;
            shape(4,:,:,3) = O;

            shape(5,:,:,1) = u;
            shape(5,:,:,2) = O;
            shape(5,:,:,3) = -s;

            shape(6,:,:,1) = O;
            shape(6,:,:,2) = u;
            shape(6,:,:,3) = -t;
        end

        function deriv = evaluateShapeDerivatives(obj, xV)
            ngaus = size(xV,2);
            nelem = size(xV,3);
            deriv = zeros(obj.ndime,obj.nnode,ngaus,nelem);
            s = xV(1,:,:);
            I = ones(size(s));
            O = zeros(size(s));

            deriv(1,1,:,:,1) = O;
            deriv(2,1,:,:,1) = I;
            deriv(3,1,:,:,1) = I;
            deriv(1,2,:,:,1) = O;
            deriv(2,2,:,:,1) = -I;
            deriv(3,2,:,:,1) = O;
            deriv(1,3,:,:,1) = O;
            deriv(2,3,:,:,1) = O;
            deriv(3,3,:,:,1) = -I;
            deriv(1,4,:,:,1) = O;
            deriv(2,4,:,:,1) = -I;
            deriv(3,4,:,:,1) = O;
            deriv(1,5,:,:,1) = O;
            deriv(2,5,:,:,1) = O;
            deriv(3,5,:,:,1) = -I;
            deriv(1,6,:,:,1) = O;
            deriv(2,6,:,:,1) = O;
            deriv(3,6,:,:,1) = O;

            deriv(1,1,:,:,2) = -I;
            deriv(2,1,:,:,2) = O;
            deriv(3,1,:,:,2) = O;
            deriv(1,2,:,:,2) = I;
            deriv(2,2,:,:,2) = O;
            deriv(3,2,:,:,2) = I;
            deriv(1,3,:,:,2) = O;
            deriv(2,3,:,:,2) = O;
            deriv(3,3,:,:,2) = -I;
            deriv(1,4,:,:,2) = I;
            deriv(2,4,:,:,2) = O;
            deriv(3,4,:,:,2) = O;
            deriv(1,5,:,:,2) = O;
            deriv(2,5,:,:,2) = O;
            deriv(3,5,:,:,2) = O;
            deriv(1,6,:,:,2) = O;
            deriv(2,6,:,:,2) = O;
            deriv(3,6,:,:,2) = -I;

            deriv(1,1,:,:,3) = -I;
            deriv(2,1,:,:,3) = O;
            deriv(3,1,:,:,3) = O;
            deriv(1,2,:,:,3) = O;
            deriv(2,2,:,:,3) = -I;
            deriv(3,2,:,:,3) = O;
            deriv(1,3,:,:,3) = I;
            deriv(2,3,:,:,3) = I;
            deriv(3,3,:,:,3) = O;
            deriv(1,4,:,:,3) = O;
            deriv(2,4,:,:,3) = O;
            deriv(3,4,:,:,3) = O;
            deriv(1,5,:,:,3) = I;
            deriv(2,5,:,:,3) = O;
            deriv(3,5,:,:,3) = O;
            deriv(1,6,:,:,3) = O;
            deriv(2,6,:,:,3) = I;
            deriv(3,6,:,:,3) = O;
        end

    end

end
