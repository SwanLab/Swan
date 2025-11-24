classdef QuadrilateralQ8 < Interpolation

    methods (Access = public)

        function obj = QuadrilateralQ8(cParams)
            obj.init(cParams);
        end

    end

    methods (Access = protected)

        function computeParams(obj)
            obj.type   = 'QUADRILATERAL_Q8';
            obj.ndime  = 2;
            obj.nnode  = 8;
            obj.pos_nodes = [...
                -1,-1;  % 1
                 1,-1;  % 2
                 1, 1;  % 3
                -1, 1;  % 4
                 0,-1;  % 5
                 1, 0;  % 6
                 0, 1;  % 7
                -1, 0]; % 8
        end

        function shape = evaluateShapeFunctions(obj,xV)
            ngaus = size(xV,2);
            nelem = size(xV,3);
            shape = zeros(obj.nnode, ngaus, nelem);
            for igaus = 1:ngaus
                s = xV(1,igaus);
                t = xV(2,igaus);
                N1 = 0.25*(1-s)*(1-t)*(-s-t-1);
                N2 = 0.25*(1+s)*(1-t)*( s-t-1);
                N3 = 0.25*(1+s)*(1+t)*( s+t-1);
                N4 = 0.25*(1-s)*(1+t)*(-s+t-1);
                N5 = 0.5*(1-s^2)*(1-t);
                N6 = 0.5*(1+s)*(1-t^2);
                N7 = 0.5*(1-s^2)*(1+t);
                N8 = 0.5*(1-s)*(1-t^2);
                shape(:,igaus,:) = [N1; N2; N3; N4; N5; N6; N7; N8];
            end
        end

        function deriv = evaluateShapeDerivatives(obj,xV)
            ngaus = size(xV,2);
            nelem = size(xV,3);
            deriv = zeros(obj.ndime, obj.nnode, ngaus, nelem);
            for igaus = 1:ngaus
                s = xV(1,igaus);
                t = xV(2,igaus);
                
                % dN/ds
                dNds = [...
                    0.25*((t-1)*(2*s+t));             % N1
                    0.25*((1-t)*(2*s-t));             % N2
                    0.25*((1+t)*(2*s+t));             % N3
                    0.25*((-1-t)*(2*s-t));            % N4
                    -s*(1-t);                         % N5
                    0.5*(1-t^2);                      % N6
                    -s*(1+t);                         % N7
                    -0.5*(1-t^2)];                    % N8
                
                % dN/dt
                dNdt = [...
                    0.25*((s-1)*(s+2*t));             % N1
                    0.25*((1+s)*(s-2*t));             % N2
                    0.25*((1+s)*(s+2*t));             % N3
                    0.25*((1-s)*(2*t-s));             % N4
                    -0.5*(1-s^2);                     % N5
                    -(1+s)*t;                         % N6
                    0.5*(1-s^2);                      % N7
                    -(1-s)*t];                        % N8

                deriv(:,:,igaus,:) = [dNds'; dNdt'];
            end
        end

    end

end
