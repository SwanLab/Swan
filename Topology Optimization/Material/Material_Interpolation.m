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
        function obj=create(ConstitutiveProperties, TypeOfMaterial, interpolation, pdim)
            switch TypeOfMaterial
                case 'ISOTROPIC'
                    switch interpolation
                        case 'SIMPALL'
                            switch pdim
                                case '2D'
                                    obj=Material_Interpolation_ISO_SIMPALL_2D(ConstitutiveProperties);
                                case '3D'
                                    obj=Material_Interpolation_ISO_SIMPALL_3D(ConstitutiveProperties);
                            end
                        case 'SIMP_Adaptative'
                            obj=Material_Interpolation_ISO_SIMP_Adaptative(ConstitutiveProperties);
                        case 'SIMP_P3'
                            obj=Material_Interpolation_ISO_SIMP_P3(ConstitutiveProperties);
                        otherwise
                            disp('Method not added')
                    end
            end
        end
        computeMatProp()
    end
end