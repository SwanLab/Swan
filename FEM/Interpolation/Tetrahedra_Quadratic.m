classdef Tetrahedra_Quadratic < Interpolation

    methods (Access = public)

        function obj = Tetrahedra_Quadratic(cParams)
            obj.init(cParams);
        end

    end

    methods (Access = protected)

        function computeParams(obj)
            obj.ndime = 3;
            obj.nnode = 10;
            obj.pos_nodes = [0 0 0;
                1 0 0;
                0 1 0;
                0 0 1;
                0.5 0 0;
                0 0.5 0;
                0 0 0.5;
                0.5 0.5 0;
                0.5 0 0.5;
                0 0.5 0.5
                ];
            % obj.isoDv = 1/6;
        end

        function shape = evaluateShapeFunctions(obj,xV)
            ngaus = size(xV,2);
            nelem = size(xV,3);
            s = xV(1,:,:);
            t = xV(2,:,:);
            u = xV(3,:,:);
            shape = zeros(obj.nnode,ngaus,nelem);
            shape(1,:,:)  = 2*s.^2 + 4.*s.*t + 4.*s.*u - 3.*s + 2.*t.^2 + 4.*t.*u - 3.*t + 2.*u.^2 -3.*u +1;
            shape(2,:,:)  = 2.*s.^2 - s;
            shape(3,:,:)  = 2.*t.^2 - t;
            shape(4,:,:)  = 2.*u.^2 - u;
            shape(5,:,:)  = 4.*s.*(-s - t - u + 1);
            shape(6,:,:)  = 4.*t.*(-s - t - u + 1);
            shape(7,:,:)  = 4.*u.*(-s - t - u + 1);
            shape(8,:,:)  = 4.*s.*t;
            shape(9,:,:)  = 4.*s.*u;
            shape(10,:,:) = 4.*t.*u;
        end

        function deriv = evaluateShapeDerivatives(obj,xV)
            ngaus = size(xV,2);
            nelem = size(xV,3);
            s = xV(1,:,:);
            t = xV(2,:,:);
            u = xV(3,:,:);
            deriv = zeros(obj.ndime,obj.nnode,ngaus,nelem);
            deriv(1,1,:,:)  = 4*s + 4*t + 4*u - 3;
            deriv(1,2,:,:)  = 4*s - 1;
            deriv(1,3,:,:)  = 0;
            deriv(1,4,:,:)  = 0;
            deriv(1,5,:,:)  = 4 - 4*t - 4*u - 8*s;
            deriv(1,6,:,:)  = -4*t;
            deriv(1,7,:,:)  = -4*u;
            deriv(1,8,:,:)  = 4*t;
            deriv(1,9,:,:)  = 4*u;
            deriv(1,10,:,:) = 0;

            deriv(2,1,:,:)  = 4*s + 4*t + 4*u - 3;
            deriv(2,2,:,:)  = 0;
            deriv(2,3,:,:)  = 4*t - 1;
            deriv(2,4,:,:)  = 0;
            deriv(2,5,:,:)  = -4*s;
            deriv(2,6,:,:)  = 4 - 8*t - 4*u - 4*s;
            deriv(2,7,:,:)  = -4*u;
            deriv(2,8,:,:)  = 4*s;
            deriv(2,9,:,:)  = 0;
            deriv(2,10,:,:) = 4*u;

            deriv(3,1,:,:)  = 4*s + 4*t + 4*u - 3;
            deriv(3,2,:,:)  = 0;
            deriv(3,3,:,:)  = 0;
            deriv(3,4,:,:)  = 4*u - 1;
            deriv(3,5,:,:)  = -4*s;
            deriv(3,6,:,:)  = -4*t;
            deriv(3,7,:,:)  = 4 - 4*t - 8*u - 4*s;
            deriv(3,8,:,:)  = 0;
            deriv(3,9,:,:)  = 4*s;
            deriv(3,10,:,:) = 4*t;
        end

    end

end
