classdef Projector < handle

    properties (Access = protected)

    end

    methods (Static, Access = public)

        function obj = create(cParams)
            p = ProjectorFactory();
            obj = p.create(cParams);
        end

    end

    methods (Static, Access = protected)

        function ord = determineQuadratureOrder(fun)
            switch fun.fType
                case 'L2'
                    ord = 'QUADRATIC';
                case 'FE'
                    ord = 'LINEAR';
                case 'GAUSSPOINTS'
                    ord = fun.quadrature.order;
            end
        end

    end

end