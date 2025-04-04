classdef TriangleQuadratic < Interpolation

    methods (Access = public)

        function obj = TriangleQuadratic(cParams)
            obj.init(cParams);
        end

    end

    methods (Access = protected)

        function computeParams(obj)
            obj.type = 'TRIANGLE_QUADRATIC';
            obj.ndime = 2;
            obj.nnode = 6;
            obj.pos_nodes = [0,0 ; 1 0; 0,1 ; 0.5,0 ; 0.5,0.5 ; 0,0.5];
        end

        function shape = evaluateShapeFunctions(obj,xV)
            ngaus = size(xV,2);
            nelem = size(xV,3);
            shape = zeros(obj.nnode, ngaus, nelem);
            for igaus=1:ngaus
                s = xV(1,igaus);
                t = xV(2,igaus);
                shape(:,igaus,:) = [(1.0-s-t)*(1.0-2*s-2*t);...
                    s*(2*s-1.0);...
                    t*(2*t-1.0);...
                    4*s*(1.0-s-t);...
                    4*s*t;...
                    4*t*(1.0-s-t)];
            end
        end

        function deriv = evaluateShapeDerivatives(obj,xV)
            ngaus = size(xV,2);
            nelem = size(xV,3);
            deriv = zeros(obj.ndime, obj.nnode, ngaus, nelem);
            for igaus=1:ngaus
                s = xV(1,igaus);
                t = xV(2,igaus);
                deriv(:,:,igaus,:) = [4*(s+t)-3,    4*s-1 , 0.0     4*(1.0-t)-8*s ,  4*t ,    -4*t; ...
                    4*(s+t)-3.0,  0.0  ,  4*t-1.0 ,    -4*s    ,   4*s  ,   4*(1.0-s)-8*t];
            end
        end

    end
end