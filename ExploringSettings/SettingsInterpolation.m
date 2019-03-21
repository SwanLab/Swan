classdef SettingsInterpolation < AbstractSettings
   
    properties (Access = public)
        constitutiveProperties = struct
        typeOfMaterial = 'ISOTROPIC'
        interpolation = 'SIMPALL'
        dim = '2D'
    end
    
    methods (Access = public)
        function obj = SettingsInterpolation()            
            obj.constitutiveProperties.rho_plus  = 1.0;
            obj.constitutiveProperties.rho_minus = 0.0;
            obj.constitutiveProperties.E_plus    = 1.0;
            obj.constitutiveProperties.E_minus   = 1e-3;
            obj.constitutiveProperties.nu_plus   = 1/3;
            obj.constitutiveProperties.nu_minus  = 1/3;
        end
    end
end