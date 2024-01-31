classdef ProjectorFactory < handle

    methods (Static, Access = public)

        function obj = create(cParams)
            switch cParams.projectorType
                case 'P0'
                    obj = Projector_toP0(cParams);
                case 'P1'
                    obj = Projector_toP1(cParams);
                case 'P2'
                    obj = Projector_toP2(cParams);
                case 'P1D'
                    obj = Projector_toP1Discontinuous(cParams);
            end
        end

    end
end