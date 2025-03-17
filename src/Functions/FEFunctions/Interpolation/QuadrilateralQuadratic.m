classdef QuadrilateralQuadratic < Interpolation

    methods (Access = public)

        function obj = QuadrilateralQuadratic(cParams)
            obj.init(cParams);
        end

    end

    methods (Access = protected)

        function computeParams(obj)
            obj.type = 'QUADRILATERAL_QUADRATIC';
            obj.ndime = 2;
            obj.nnode = 9;
            obj.pos_nodes = [-1,-1 ; 1 -1 ; 1,1 ; -1,1 ; 0,-1 ; 1,0 ; 0,1 ; -1,0 ; 0,0];
        end

        function shape = evaluateShapeFunctions(obj,xV)
            ngaus = size(xV,2);
            nelem = size(xV,3);
            shape = zeros(obj.nnode, ngaus, nelem);
            for igaus=1:ngaus
                s = xV(1,igaus);
                t = xV(2,igaus);
                shape(:,igaus,:) = [(s^2*t^2)/4 - (s^2*t)/4 - (s*t^2)/4 + (s*t)/4;
                    (s^2*t^2)/4 - (s^2*t)/4 + (s*t^2)/4 - (s*t)/4;
                    (s^2*t^2)/4 + (s^2*t)/4 + (s*t^2)/4 + (s*t)/4;
                    (s^2*t^2)/4 + (s^2*t)/4 - (s*t^2)/4 - (s*t)/4;
                    - (s^2*t^2)/2 + (s^2*t)/2 + t^2/2 - t/2;
                    s/2 - (s^2*t^2)/2 - (s*t^2)/2 + s^2/2;
                    t/2 - (s^2*t^2)/2 - (s^2*t)/2 + t^2/2;
                    - (s^2*t^2)/2 + s^2/2 + (s*t^2)/2 - s/2;
                    s^2*t^2 - s^2 - t^2 + 1];
            end
        end

        function deriv = evaluateShapeDerivatives(obj,xV)
            ngaus = size(xV,2);
            nelem = size(xV,3);
            deriv = zeros(obj.ndime, obj.nnode, ngaus, nelem);
            for igaus=1:ngaus
                s = xV(1,igaus);
                t = xV(2,igaus);
                deriv(:,:,igaus,:) = [
                    t/4 - (s*t)/2 + (s*t^2)/2 - t^2/4,...
                    (s*t^2)/2 - (s*t)/2 - t/4 + t^2/4,...
                    t/4 + (s*t)/2 + (s*t^2)/2 + t^2/4,...
                    (s*t)/2 - t/4 + (s*t^2)/2 - t^2/4,...
                    - s*t^2 + s*t,...
                    s - s*t^2 - t^2/2 + 1/2,...
                    - s*t^2 - s*t,...
                    s - s*t^2 + t^2/2 - 1/2,...
                    2*s*t^2 - 2*s;

                    s/4 - (s*t)/2 + (s^2*t)/2 - s^2/4,...
                    (s*t)/2 - s/4 + (s^2*t)/2 - s^2/4,...
                    s/4 + (s*t)/2 + (s^2*t)/2 + s^2/4,...
                    (s^2*t)/2 - (s*t)/2 - s/4 + s^2/4,...
                    t - s^2*t + s^2/2 - 1/2,...
                    - t*s^2 - t*s,...
                    t - s^2*t - s^2/2 + 1/2,...
                    - t*s^2 + t*s,...
                    2*t*s^2 - 2*t];
            end
        end

    end

end