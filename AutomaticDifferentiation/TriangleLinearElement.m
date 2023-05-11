classdef TriangleLinearElement

    properties

        coordElem

    end

    methods (Access = public)

        function obj = TriangleLinearElement(a) %triangleElement class constructor;

            obj.coordElem = a; %parametric coordinates

        end


        function [J,N] = assembleJ(obj) %Assembly Jacobian Matrix.
            %coordinates of the triangle vertex (x1, y1, x2, y2, etc) (random)coordElem

            coef = zeros(2,3);

            % matrix of coefficients
            for i = 1 : 2 %coefficients of Nj *  parametric coordinates. Linear equation

                coef(i,1) = obj.coordElem(2,i) - obj.coordElem(1,i);
                coef(i,2) = obj.coordElem(3,i) - obj.coordElem(1,i);
                coef(i,3) = obj.coordElem(1,i);

            end

            %AD using my code
            xi = ValDerForward(0,[1 0]);
            eta = ValDerForward(0,[0 1]);

            coord(1) = coef(1,1) * xi + coef(1,2) * eta + coef(1,3);
            coord(2) = coef(2,1) * xi + coef(2,2) * eta + coef(2,3);

            N = [ coord(1).double; coord(2).double ];

            %Jacobian
            J = N(:,2:end);

        end
    end
end
