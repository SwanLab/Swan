classdef FilterFactory < handle

    methods (Access = public, Static)

        function filter = create(cParams)
            switch cParams.filterType
                case 'P1'
                    switch cParams.designVarType
                        case {'Density','MicroParams'}
                            filter = Filter_P1_Density(cParams);
                        case 'LevelSet'
                            filter = Filter_P1_LevelSet(cParams);
                    end
                case 'PDE'
                    switch cParams.designVarType
                        case {'Density','MicroParams'}
                            filter = Filter_PDE_Density(cParams);
                        case 'LevelSet'
                            filter = Filter_PDE_LevelSet(cParams);
                    end
            end
        end

    end

end