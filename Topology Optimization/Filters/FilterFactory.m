classdef FilterFactory < handle

    % Refactoring ToDo:
% - Test P0Function in unfitted not available
% - 'evaluate' method of CharacteristicFunction inside RHSUnfitted?

    methods (Access = public, Static)

        function filter = create(cParams)
            switch cParams.filterType
                case 'P1'
                    switch cParams.designVarType
                        case {'Continuous','Density','MicroParams'}
                            filter = FilterP1(cParams);
                        case 'LevelSet'
                            filter = FilterP1Unfitted(cParams);
                    end
                case 'PDE'
                    switch cParams.designVarType
                        case {'Continuous','Density','MicroParams'}
                            filter = FilterPDE(cParams);
                        case 'LevelSet'
                            filter = FilterPDEUnfitted(cParams);
                    end
            end
        end

    end

end