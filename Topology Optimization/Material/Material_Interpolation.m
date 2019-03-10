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
        function obj=create(cParams)
            switch cParams.TypeOfMaterial
                case 'ISOTROPIC'
                    switch cParams.Interpolation
                        case 'SIMPALL'
                            switch cParams.Dim
                                case '2D'
                                    obj=Material_Interpolation_ISO_SIMPALL_2D(cParams.ConstitutiveProperties);
                                case '3D'
                                    obj=Material_Interpolation_ISO_SIMPALL_3D(cParams.ConstitutiveProperties);
                            end
                        case 'SIMP_Adaptative'
                            obj=Material_Interpolation_ISO_SIMP_Adaptative(cParams.ConstitutiveProperties);
                        case 'SIMP_P3'
                            obj=Material_Interpolation_ISO_SIMP_P3(cParams.ConstitutiveProperties);
                        otherwise
                            disp('Method not added')
                    end
            end
        end
        computeMatProp()
    end
end