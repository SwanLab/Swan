classdef SettingsConstraint < SettingsCC
    
    properties (Access = protected)
        defaultParamsName = 'paramsConstraint.json'
    end
    
    properties (Access = public)
        dualVariable
    end
    
    methods (Access = public)
        
        function obj = SettingsConstraint(varargin)
            if nargin == 1
                obj.loadParams(varargin{1});
            end
        end
        
    end
    
    methods (Access = public, Static)
        
        function s = create(cParams,settings)

        end
        
    end
    
end