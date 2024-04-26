classdef ProjectorFactory < handle

    methods (Static, Access = public)

        function obj = create(cParams)
            switch cParams.projectorType
                case 'P0'
                    obj = Projector_toP0(cParams);
                case {'P1','P2','P3'}
                    obj = Projector_toLagrangian(cParams);
                case 'P1D'
                    obj = Projector_toP1Discontinuous(cParams);
            end
        end

    end
end