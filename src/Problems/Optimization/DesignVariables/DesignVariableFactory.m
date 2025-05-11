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
                case 'DensityAndBound'
                    designVar = DensityAndBound(cParams);
                case 'MultiLevelSet'
                    designVar = MultiLevelSet(cParams);
            end
        end

    end

end

