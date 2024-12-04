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
            switch class(fun)
                case 'L2Function'
                    ord = 2;
                case 'FEFunction'
                    ord = 2; % needed to project P1 to P1D             
                otherwise
                    ord = 3;
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