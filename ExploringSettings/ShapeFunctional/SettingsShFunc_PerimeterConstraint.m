classdef SettingsShFunc_PerimeterConstraint < SettingsShapeFunctional

    properties
        Perimeter_target
    end
    
     methods (Access = public)
        
        function obj = SettingsShFunc_PerimeterConstraint(varargin)
            obj.defaultParamsName = 'paramsShapeFunctional_PerimeterConstraint';
            switch nargin
                case 0
                    obj.loadParams(obj.defaultParamsName);
                case 1
                    obj.loadParams(varargin{1});
            end
        end
        
    end
    
end