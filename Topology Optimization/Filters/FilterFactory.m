classdef FilterFactory < handle
    
    methods (Access = public)
        
        function filter = create(obj,factoryParams)
            switch factoryParams.type
                case 'P1'
                    switch factoryParams.designVar
                        case 'DENSITY'
                            filter = Filter_P1_Density(cParams);
                        case 'LEVELSET'
                            filter = Filter_P1_LevelSet(cParams);
                    end
                case 'PDE'
                    switch factoryParams.designVar
                        case 'DENSITY'
                            filter = Filter_PDE_Density(cParams);
                        case 'LEVELSET'
                            filter = Filter_PDE_LevelSet(cParams);
                    end
            end
        end
        
    end
    
end