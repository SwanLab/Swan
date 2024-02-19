classdef Quadrilateral_Quadratic < Interpolation

    methods (Access = public)

        function obj = Quadrilateral_Quadratic(cParams)
            obj.init(cParams);
            obj.computeParams();
        end

        function shape = computeShapeFunctions(obj,posgp)
            ngaus = size(posgp,2);
            shape = zeros(obj.nnode, ngaus);
            for igaus=1:ngaus
                s = posgp(1,igaus);
                t = posgp(2,igaus);
                shape(:,igaus) = [(s^2*t^2)/4 - (s^2*t)/4 - (s*t^2)/4 + (s*t)/4;
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

        function deriv = computeShapeDerivatives(obj,posgp)
            ngaus = size(posgp,2);
            nelem = size(posgp,3);
            for igaus=1:ngaus
                s = posgp(1,igaus);
                t = posgp(2,igaus);
                obj.deriv(:,:,igaus) = [
                    t/4 - (s*t)/2 + (s*t^2)/2 - t^2/4,
                    (s*t^2)/2 - (s*t)/2 - t/4 + t^2/4,
                    t/4 + (s*t)/2 + (s*t^2)/2 + t^2/4,
                    (s*t)/2 - t/4 + (s*t^2)/2 - t^2/4,
                    - s*t^2 + s*t,
                    s - s*t^2 - t^2/2 + 1/2,
                    - s*t^2 - s*t,
                    s - s*t^2 + t^2/2 - 1/2,
                    2*s*t^2 - 2*s;

                    s/4 - (s*t)/2 + (s^2*t)/2 - s^2/4,
                    (s*t)/2 - s/4 + (s^2*t)/2 - s^2/4,
                    s/4 + (s*t)/2 + (s^2*t)/2 + s^2/4,
                    (s^2*t)/2 - (s*t)/2 - s/4 + s^2/4,
                    t - s^2*t + s^2/2 - 1/2,
                    - t*s^2 - t*s,
                    t - s^2*t - s^2/2 + 1/2,
                    - t*s^2 + t*s,
                    2*t*s^2 - 2*t];
            end
        end

    end

    methods (Access = private)

        function computeParams(obj)
            obj.type = 'QUADRILATERAL_QUADRATIC';
            obj.ndime = 2;
            obj.nnode = 9;
            obj.pos_nodes = [-1,-1 ; 1 -1 ; 1,1 ; -1,1 ; 0,-1 ; 1,0 ; 0,1 ; -1,0 ; 0,0];
        end

    end

end