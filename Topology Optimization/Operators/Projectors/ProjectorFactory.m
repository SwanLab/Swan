classdef ProjectorFactory < handle

    methods (Static, Access = public)

        function obj = create(cParams)
            switch cParams.projectorType
                case 'LINEAR'
                    obj = Projector_toLagrangian(cParams);
                case 'QUADRATIC'
                    obj = Projector_toLagrangian(cParams);
                case 'CUBIC'
                    obj = Projector_toLagrangian(cParams);
                case 'ORDER4'
                    obj = Projector_toLagrangian(cParams);
                case 'P1D'
                    obj = Projector_toP1Discontinuous(cParams);
                case 'H1P1'
                    obj = H1Projector_toP1(cParams);
                case 'H1P1D'
                    obj = H1Projector_toP1Discontinuous(cParams);
            end
        end

    end
end