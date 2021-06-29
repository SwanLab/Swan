classdef FilterFactory < handle
    
    methods (Access = public)
        
        function filter = create(obj,factoryParams)
            cParams = factoryParams;
            switch factoryParams.filterType
                case 'P1'
                    switch factoryParams.designVarType                                                
                        case {'Density','MicroParams'}
                            filter = Filter_P1_Density(cParams);
                        case 'LevelSet'
                            filter = Filter_P1_LevelSet(cParams);
                    end
                case 'PDE'
                    switch factoryParams.designVarType
                        case {'Density','MicroParams'}
                            filter = Filter_PDE_Density(cParams);
                        case 'LevelSet'
                            filter = Filter_PDE_LevelSet(cParams);
                    end
            end
        end
        
    end
    
end