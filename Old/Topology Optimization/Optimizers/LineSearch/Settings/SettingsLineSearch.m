classdef SettingsLineSearch < AbstractSettings
    
    properties (Access = protected)
        defaultParamsName = 'paramsLineSearch.json'
    end
    
    properties (Access = public)
        type
        optimizerType
        filename               
        scalarProduct
        epsilon
        HJiter0
        lineSearchInitiatorSettings
        designVariable
        objectiveFunction
        rate
        incrementFactor
    end
    
    methods (Access = public)
        
        function obj = SettingsLineSearch(varargin)
            if nargin == 1
                obj.loadParams(varargin{1});
            end
        end
        
    end
    
end