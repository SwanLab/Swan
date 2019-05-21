classdef SettingsConstraint_OLD < SettingsCC_OLD
    
    properties (Access = protected)
        defaultParamsName = 'paramsConstraint.json'
    end
    
    properties (Access = public)
        dualVariable
    end
    
    methods (Access = public)
        
        function obj = SettingsConstraint_OLD(varargin)
            if nargin == 1
                obj.loadParams(varargin{1});
            end
            obj.init();
        end
        
    end
    
end