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
                case 'SquareColumn'
                    designVar = SquareColumn(cParams);
                case 'RectangularColumn'
                    designVar = RectangularColumn(cParams);
                case 'RectangularHoleColumn'
                    designVar = RectangularHoleColumn(cParams);
                case 'HoleColumn'
                    designVar = HoleColumn(cParams);
                case 'LshapeColumn'
                    designVar = LshapeColumn(cParams);                    
                case 'LevelSetEigModes'
                    designVar = LevelSetEigModes(cParams);
            end
        end
        
    end

end

