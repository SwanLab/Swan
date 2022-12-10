classdef DesignVariableFactory < handle

    methods (Access = public, Static)
        
        function designVar = create(cParams)
            switch cParams.type
                case 'LevelSet'
                    designVar = LevelSet(cParams);
                case 'Density'
                    designVar = Density(cParams);
                case 'DensityEigModes'
                    designVar = DensityEigModes(cParams);
                case 'MicroParams'
                    designVar = MicroParams(cParams);
                case 'AreaColumn'
                    designVar = AreaColumn(cParams);
                case 'RadiusColumn'
                    designVar = RadiusColumn(cParams);
                case 'LevelSetEigModes'
                    designVar = LevelSetEigModes(cParams);
            end
        end
        
    end

end

