classdef SettingsLevelSetCreator < AbstractSettings
    
    properties (Access = protected)
        defaultParamsName = 'paramsLevelSetCreator.json'
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
        
        function obj = create(obj,s)
           factory = SettingsLevelSetCreatorFactory();
           obj     = factory.create(s);
        end
        
    end
    
end