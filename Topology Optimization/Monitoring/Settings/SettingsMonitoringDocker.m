classdef SettingsMonitoringDocker < AbstractSettings
    
    properties (Access = protected)
        defaultParamsName = 'paramsMonitoringDocker.json'
    end
    
    properties (Access = public)
        showOptParams
        refreshInterval
        shallDisplayDesignVar
        shallShowBoundaryConditions
        
        problemID
        costFuncNames
        costWeights
        constraintFuncs
        dim
        
        designVariable
        dualVariable
        optimizerName
        cost
        constraint
        
        convergenceVars
    end
    
    methods (Access = public)
        
        function obj = SettingsMonitoringDocker(varargin)
            if nargin == 1
                obj.loadParams(varargin{1});
            end
        end
        
    end
    
end