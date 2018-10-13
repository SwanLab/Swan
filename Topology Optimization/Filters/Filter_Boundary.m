classdef Filter_Boundary < Filter
    methods (Static)
        function obj = create(settings)
            switch settings.filter
                case 'P1'
                    switch settings.pdim
                        case '2D'
                            obj = Filter_P1_LevelSet_2D(settings.filename,settings.ptype);
                        case '3D'
                            obj = Filter_P1_LevelSet_3D(settings.filename,settings.ptype);
                    end
                case 'PDE'
                    switch settings.pdim
                        case '2D'
                            obj = Filter_PDE_LevelSet_2D(settings.filename,settings.ptype);
                        case '3D'
                            obj = Filter_PDE_LevelSet_3D(settings.filename,settings.ptype);
                    end
            end
        end
    end
end