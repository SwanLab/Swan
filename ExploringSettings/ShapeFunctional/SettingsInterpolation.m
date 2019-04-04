classdef SettingsInterpolation < AbstractSettings
    
    
    properties (Access = protected)
        defaultParamsName = 'paramsMaterialInterpolation'
    end
   
    properties (Access = public)
        constitutiveProperties = struct
        typeOfMaterial
        interpolation 
        dim 
    end
    
    methods (Access = public)
        function obj = SettingsInterpolation(varargin)            
            obj.constitutiveProperties.rho_plus  = 1.0;
            obj.constitutiveProperties.rho_minus = 0.0;
            obj.constitutiveProperties.E_plus    = 1.0;
            obj.constitutiveProperties.E_minus   = 1e-3;
            obj.constitutiveProperties.nu_plus   = 1/3;
            obj.constitutiveProperties.nu_minus  = 1/3;
            if nargin == 1
                obj.loadParams(varargin{1})
            end
        end
    end
end