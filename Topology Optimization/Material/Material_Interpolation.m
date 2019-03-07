classdef Material_Interpolation < handle
    properties
        E_plus
        E_minus
        nu_plus
        nu_minus
        rho_plus
        rho_minus
    end
    
    methods (Static)
        function obj=create(settings)
            switch settings.TypeOfMaterial
                case 'ISOTROPIC'
                    switch settings.Interpolation
                        case 'SIMPALL'
                            switch settings.Dim
                                case '2D'
                                    obj=Material_Interpolation_ISO_SIMPALL_2D(settings.ConstitutiveProperties);
                                case '3D'
                                    obj=Material_Interpolation_ISO_SIMPALL_3D(settings.ConstitutiveProperties);
                            end
                        case 'SIMP_Adaptative'
                            obj=Material_Interpolation_ISO_SIMP_Adaptative(settings.ConstitutiveProperties);
                        case 'SIMP_P3'
                            obj=Material_Interpolation_ISO_SIMP_P3(settings.ConstitutiveProperties);
                        otherwise
                            disp('Method not added')
                    end
            end
        end
        computeMatProp()
    end
end