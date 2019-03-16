classdef SettingsNumericalHomogenizer < DefaultSettings
  
    properties (Access = public)
        interpParams = SettingsInterpolation()
        levelSetCreatorParams = SettingsLevelSetCreator()
        volumeShFuncParams = SettingsShapeFunctional()
        elementDensityCreatorType = 'ElementalDensityCreatorByLevelSetCreator'
        outFileName = 'RVE_Square_Triangle_Fine'
        testName = 'RVE_Square_Triangle_Fine.m'
        print = false
        iter = 0
        pdim = '2D'
    end
    
    methods
        
        function obj = SettingsNumericalHomogenizer()
            obj.volumeShFuncParams.filename = 'RVE_Square_Triangle_Fine.m';
            obj.volumeShFuncParams.ptype = 'MICRO';
        end
        
    end
    
end