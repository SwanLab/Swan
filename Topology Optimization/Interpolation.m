classdef Interpolation < handle
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
                            obj=Interpolation_ISO_SIMPALL(TOL);
                        otherwise
                            disp('Method not added')
                    end
           end
        end
    computeMatProp()
    end
end