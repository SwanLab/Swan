classdef SettingsNumericalHomogenizer < DefaultSettings
  
    properties (Access = public)
        interpParams = SettingsInterpolation()
        levelSetCreatorParams = SettingsLevelSetCreator()
        volumeShFuncParams = SettingsShapeFunctional()
        elementDensityCreatorType = 'ElementalDensityCreatorByLevelSetCreator'
        outFileName = 'output'
        testName = 'please specify filename'
        print = false
        iter = 0
        pdim = '2D'
    end
    
    methods
    end
end