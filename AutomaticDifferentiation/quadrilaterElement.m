classdef QuadrilaterElement

    properties

        coordElem

    end

    methods (Access = public)

        function obj = QuadrilaterElement(a) %quadrilaterElement class constructor;

            obj.coordElem = a; %parametric coordinates

        end


        function [J,N] = assembleJ(obj) %Assembly Jacobian Matrix.

            coef = zeros(2,4);

            % matrix of coefficients
            for i = 1 : 2 %coefficients of Nj *  parametric coordinates. Linear equation

                coef(i,1) = 0.25 * ( obj.coordElem(1,i) - obj.coordElem(2,i) + obj.coordElem(3,i) - obj.coordElem(4,i) );
                coef(i,2) = 0.25 * ( - obj.coordElem(1,i) + obj.coordElem(2,i) + obj.coordElem(3,i) - obj.coordElem(4,i) );
                coef(i,3) = 0.25 * ( - obj.coordElem(1,i) - obj.coordElem(2,i) + obj.coordElem(3,i) + obj.coordElem(4,i) );
                coef(i,4) = 0.25 * ( obj.coordElem(1,i) + obj.coordElem(2,i) + obj.coordElem(3,i) + obj.coordElem(4,i) );

            end

            %AD using my code
            xi = ValDerForward(1,[1 0]);
            eta = ValDerForward(1,[0 1]);

            coord(1) = coef(1,1) * xi * eta + coef(1,2) * xi + coef(1,3) * eta + coef(1,4);
            coord(2) = coef(2,1) * xi * eta + coef(2,2) * xi + coef(2,3) * eta + coef(2,4);

            N = [ coord(1).double; coord(2).double ];

            %Jacobian
            J = N(:,2:end);
        end
    end
end