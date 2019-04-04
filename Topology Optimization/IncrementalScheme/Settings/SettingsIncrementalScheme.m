classdef SettingsIncrementalScheme < AbstractSettings
    
    properties (Access = protected)
        defaultParamsName = 'paramsIncrementalScheme'
    end
    
    properties (Access = public)
        nSteps
        settingsTargetParams        
        shallPrintIncremental        
        mesh
    end
    
    methods (Access = public)
        
        function obj = SettingsIncrementalScheme(varargin)
            if nargin == 1
                obj.loadParams(varargin{1});
            end
        end
        
    end
    
end