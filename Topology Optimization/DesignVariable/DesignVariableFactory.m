classdef DesignVariableFactory < handle

    methods (Access = public, Static)
        
        function designVar = create(cParams)
            switch cParams.type
                case 'LevelSet'
                    designVar = LevelSet(cParams);
                case 'Density'
                    designVar = Density(cParams);
                case 'MicroParams'
                    designVar = MicroParams(cParams);
                case 'AreaColumn'
                    designVar = AreaColumn(cParams);
                case 'TrussDiscrete'
                    designVar = TrussDiscrete(cParams);
            end
        end
        
    end

end

