classdef Triangle_Quadratic < Interpolation

    methods (Access = public)

        function obj = Triangle_Quadratic(cParams)
            obj.init(cParams);
            obj.computeParams();
        end

        function shape = computeShapeFunctions(obj,posgp)
            ngaus = size(posgp,2);
            shape = zeros(obj.nnode, ngaus);
            for igaus=1:ngaus
                s = posgp(1,igaus);
                t = posgp(2,igaus);
                shape(:,igaus) = [(1.0-s-t)*(1.0-2*s-2*t);...
                    s*(2*s-1.0);...
                    t*(2*t-1.0);...
                    4*s*(1.0-s-t);...
                    4*s*t;...
                    4*t*(1.0-s-t)];
            end
        end

        function deriv = computeShapeDerivatives(obj,posgp)
            ngaus = size(posgp,2);
            deriv = zeros(obj.ndime, obj.nnode, ngaus);
            for igaus=1:ngaus
                s = posgp(1,igaus);
                t = posgp(2,igaus);
                deriv(:,:,igaus) = [4*(s+t)-3,    4*s-1 , 0.0     4*(1.0-t)-8*s ,  4*t ,    -4*t; ...
                    4*(s+t)-3.0,  0.0  ,  4*t-1.0 ,    -4*s    ,   4*s  ,   4*(1.0-s)-8*t];
            end
        end

    end

    methods (Access = private)

        function computeParams(obj)
            obj.type = 'TRIANGLE_QUADRATIC';
            obj.ndime = 2;
            obj.nnode = 6;
            obj.pos_nodes = [0,0 ; 1 0; 0,1 ; 0.5,0 ; 0.5,0.5 ; 0,0.5];
        end

    end
end