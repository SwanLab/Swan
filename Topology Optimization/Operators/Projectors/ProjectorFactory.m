classdef ProjectorFactory < handle

    methods (Static, Access = public)

        function obj = create(cParams)
            switch cParams.projectorType
                case 'toP0'
                    obj = Projector_toP0(cParams);
                case 'toP1'
                    obj = Projector_toP1(cParams);
                case 'toP1Disc'
                    obj = Projector_toP1Discontinuous(cParams);
                case 'toH1P1'
                    obj = H1Projector_toP1(cParams);
                case 'toH1P1Disc'
                    obj = H1Projector_toP1Discontinuous(cParams);
            end
        end

    end
end