classdef Tetrahedra_Quadratic < Interpolation

    methods (Access = public)

        function obj = Tetrahedra_Quadratic(cParams)
            obj.init(cParams);
            obj.computeParams();
            obj.computeCases();
        end

        function shape = computeShapeFunctions(obj,posgp)
            ngaus = size(posgp,2);
            nelem = size(posgp,3);
            s = posgp(1,:,:);
            t = posgp(2,:,:);
            u = posgp(3,:,:);
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

        function deriv = computeShapeDerivatives(obj,posgp)
            ngaus = size(posgp,2);
            nelem = size(posgp,3);
            s = posgp(1,:,:);
            t = posgp(2,:,:);
            u = posgp(3,:,:);
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

    methods (Access = private)

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
            obj.isoDv = 1/6;
        end

        function computeCases(obj)
            obj.iteration = [1 1 1 2 2 3;
                2 3 4 3 4 4];
            obj.cases(:,:,1) = [7 6 5 1
                5 6 7 4
                5 6 4 2
                4 3 2 6;
                zeros(2,4)];
            obj.cases(:,:,2) = [5 6 7 2
                7 6 5 3
                5 7 3 1
                7 3 1 4;
                zeros(2,4)];
            obj.cases(:,:,3) = [7 6 5 3
                5 6 7 2
                2 7 5 1
                2 7 1 4;
                zeros(2,4)];
            obj.cases(:,:,4) = [5 6 7 4
                7 6 5 3
                6 5 3 1
                6 1 3 2;
                zeros(2,4)];
            obj.cases(:,:,5) = [2 5 7 3
                5 8 7 3
                5 6 8 3
                7 8 1 4
                7 6 5 1
                7 8 6 1];
            obj.cases(:,:,6) = [6 8 3 4
                6 8 7 3
                7 5 6 3
                6 5 2 1
                8 5 7 2
                8 6 5 2];
            obj.cases(:,:,7) = [2 8 6 4
                5 8 6 2
                5 7 8 2
                6 3 5 1
                8 7 5 3
                6 8 5 3];
            obj.main_loop = [4 4];
            obj.extra_cases = [5,6,7];
            obj.selectcases = [1 0 0;
                2 0 0;
                3 6 0;
                4 7 0;
                0 5 0;
                0 7 4;
                0 6 3
                0 0 2
                0 0 1];
        end
    end

end
