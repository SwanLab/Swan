classdef DesignVariableFactory < handle

    methods (Access = public, Static)

        function designVar = create(cParams)
            switch cParams.type
                case 'LevelSet'
                    designVar = LevelSet(cParams);
                case 'Density'
                    designVar = Density(cParams);
                case 'MultiField'
                    designVar = MultiField(cParams);
                case 'AreaColumn'
                    designVar = AreaColumn(cParams);
                case 'DensityAndBound'
                    designVar = DensityAndBound(cParams);
                case 'MultiLevelSet'
                    designVar = MultiLevelSet(cParams);
                case 'Damage'
                    designVar = Damage(cParams);
                otherwise
                    error('DesignVariableFactory:UnknownType', ...
                        'Unknown design variable type: %s', cParams.type);
            end
        end

    end

end