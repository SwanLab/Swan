classdef ProjectorFactory < handle

    methods (Static, Access = public)

        function obj = create(cParams)
            switch cParams.projectorType
                case {'P0','P1','P1D','P2','P3'}
                    obj = ProjectorToLagrangian(cParams);
                case 'RT'
                    obj = Projector_toRaviartThomas(cParams);
                case 'N'
                    obj = Projector_toNedelec(cParams);
            end
        end

    end
end