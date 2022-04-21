classdef NewFilterFactory < handle

    methods (Access = public, Static)

        function filter = create(factoryParams)
            cParams = factoryParams;
            switch factoryParams.filterType
                case 'P1'
                    switch factoryParams.designVarType
                        case {'Density','MicroParams'}
                            filter = NewFilter_P1_Density(cParams);
                        case 'LevelSet'
                            filter = NewFilter_P1_LevelSet(cParams);
                    end
                case 'PDE'
                    switch factoryParams.designVarType
                        case {'Density','MicroParams'}
                            filter = NewFilter_PDE_Density(cParams);
                        case 'LevelSet'
                            filter = NewFilter_PDE_LevelSet(cParams);
                    end
            end
        end

    end

end