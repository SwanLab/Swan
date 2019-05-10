classdef SettingsLevelSetCreator < AbstractSettings
    
    properties (Access = protected)
        defaultParamsName = 'paramsLevelSetCreator'
    end
    
    properties (Access = public)
        type
        ndim
        coord
        value
    end
    
    methods (Access = public)
        
        function obj = SettingsLevelSetCreator(varargin)
            if nargin == 1
                obj.loadParams(varargin{1})
            end
        end
        
    end
    
end