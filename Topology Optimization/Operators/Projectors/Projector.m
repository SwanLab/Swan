classdef Projector < handle

    properties (Access = protected)
        mesh
        connec
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
%                     ord = 'LINEAR';
                    ord = 'QUADRATIC'; % needed to project P1 to P1D
                case 'GAUSSPOINTS'
                    ord = fun.quadrature.order;
            end
        end

    end

    methods (Access = protected)

        function init(obj, cParams)
            obj.mesh   = cParams.mesh;
            obj.connec = cParams.mesh.connec;
        end

    end

end