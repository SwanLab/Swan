classdef SettingsInterpolation < AbstractSettings
    
    properties (Access = protected)
        defaultParamsName = 'paramsMaterialInterpolation.json'
    end
    
    properties (Access = public)
        constitutiveProperties
        typeOfMaterial
        interpolation
        dim
        ndim
        type
        nElem
        ngaus
        simpAllType
    end
    
    methods (Access = public)
        
        function obj = SettingsInterpolation(varargin)
            if nargin == 1
                obj.loadParams(varargin{1})
            end
            obj.init();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            s = obj.constitutiveProperties;
            if isstruct(s)
                obj.constitutiveProperties = SettingsConstitutiveProperties(s);
            end
        end
        
    end
    
end