classdef Filter_Boundary < Filter
    methods (Static)
        function obj = create(settings)
            switch settings.filter
                case 'P1'
                    obj = Filter_P1_LevelSet();
                case 'PDE'
                    obj = Filter_PDE_LevelSet();
            end
        end
    end
end