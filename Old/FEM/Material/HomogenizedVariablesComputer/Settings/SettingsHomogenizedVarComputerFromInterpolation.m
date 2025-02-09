classdef SettingsHomogenizedVarComputerFromInterpolation < ....
        SettingsHomogenizedVarComputer
    
    properties (Access = protected)
        defaultParamsName = 'paramsHomogenizedVarComputerFromInterpolation.json'
    end
    
    properties (Access = public)
        interpolation
        typeOfMaterial
        constitutiveProperties
        nelem
        ptype
        interpolationSettings
    end
    
    methods (Access = public)
        
        function obj = SettingsHomogenizedVarComputerFromInterpolation(varargin)
            if nargin == 1
                obj.loadParams(varargin{1})
            end
            obj.init();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.initConstitutiveProperties();
            obj.initInterpolationSettings();
        end
        
        function initConstitutiveProperties(obj)
            s = obj.constitutiveProperties;
            obj.constitutiveProperties = SettingsConstitutiveProperties(s);
        end
        
        function initInterpolationSettings(obj)
            s.interpolation           = obj.interpolation;
            s.dim                     = obj.dim;
            s.typeOfMaterial          = obj.typeOfMaterial;
            s.constitutiveProperties  = obj.constitutiveProperties;
            obj.interpolationSettings = SettingsInterpolation(s);
        end
        
    end
    
end