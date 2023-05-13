classdef BoundaryConditionsFactory < handle

    methods (Access = public, Static)
        function obj = create(cParams)
            switch cParams.solMode
                case 'DISP'
                    obj = BoundaryConditionsLineal(cParams);
                case 'FLUC'
                    obj = BoundaryConditionsPeriodic(cParams);
            end
        end

    end
end