classdef SettingsCost < SettingsCC
    
    properties (Access = protected)
        defaultParamsName = 'paramsCost.json'
    end
    
    properties (Access = public)
        weights
    end
    
    methods (Access = public)
        
        function obj = SettingsCost(varargin)
            if nargin == 1
                obj.loadParams(varargin{1});
            end
            obj.init();
        end

    end
    
end