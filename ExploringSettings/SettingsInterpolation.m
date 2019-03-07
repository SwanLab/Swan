classdef SettingsInterpolation < DefaultSettings
    properties
        ConstitutiveProperties=struct
        TypeOfMaterial = 'ISOTROPIC'
        Interpolation = 'SIMPALL'
        Dim = '2D'
    end
    methods
        function obj = SettingsInterpolation()            
            obj.ConstitutiveProperties.rho_plus  = 1.0;
            obj.ConstitutiveProperties.rho_minus = 0.0;
            obj.ConstitutiveProperties.E_plus    = 1.0;
            obj.ConstitutiveProperties.E_minus   = 1e-3;
            obj.ConstitutiveProperties.nu_plus   = 1/3;
            obj.ConstitutiveProperties.nu_minus  = 1/3;
        end
    end
end