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
        function obj=create(TOL, material, method)
            switch material
                case 'ISOTROPIC'
                    switch method
                        case 'SIMPALL'
                            obj=Material_Interpolation_ISO_SIMPALL(TOL);
                        case 'SIMP_Adaptative'
                            obj=Material_Interpolation_ISO_SIMP_Adaptative(TOL);
                        case 'SIMP_P3'
                            obj=Material_Interpolation_ISO_SIMP_P3(TOL);
                        otherwise
                            disp('Method not added')
                    end
            end
        end
        computeMatProp()
    end
end