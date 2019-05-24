classdef SettingsLineSearch < AbstractSettings
    
    properties (Access = protected)
        defaultParamsName = 'paramsLineSearch.json'
    end
    
    properties (Access = public)
        type
        optimizerType
        filename
        kfrac
        kappaMultiplier
        scalarProduct
        epsilon
        HJiter0
    end
    
    methods (Access = public)
        
        function obj = SettingsLineSearch(varargin)
            if nargin == 1
                obj.loadParams(varargin{1});
            end
        end
        
    end
    
end