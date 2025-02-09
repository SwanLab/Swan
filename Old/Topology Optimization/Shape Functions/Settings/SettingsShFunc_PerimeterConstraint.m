classdef SettingsShFunc_PerimeterConstraint < SettingsShapeFunctional
    
    properties (Access = public)
        perimeterTarget
    end
    
    methods (Access = public)
        
        function obj = SettingsShFunc_PerimeterConstraint(varargin)
            obj.defaultParamsName = 'paramsShapeFunctional_PerimeterConstraint';
            obj.type = 'perimeterConstraint';
            if nargin == 1
                    obj.loadParams(varargin{1});
            end
        end
        
    end
    
end