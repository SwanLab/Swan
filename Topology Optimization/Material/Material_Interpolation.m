classdef Material_Interpolation < handle
    
    properties (SetAccess = protected, GetAccess = public)
        E_plus
        E_minus
        nu_plus
        nu_minus
        rho_plus
        rho_minus
    end
    
    methods (Static)
        function obj = create(cParams)
            constParams = cParams.constitutiveProperties;
            switch cParams.typeOfMaterial
                case 'ISOTROPIC'
                    switch cParams.interpolation
                        case 'SIMPALL'
                            switch cParams.dim
                                case '2D'
                                    obj=Material_Interpolation_ISO_SIMPALL_2D(constParams);
                                case '3D'
                                    obj=Material_Interpolation_ISO_SIMPALL_3D(constParams);
                            end
                        case 'SIMP_Adaptative'
                            obj=Material_Interpolation_ISO_SIMP_Adaptative(constParams);
                        case 'SIMP_P3'
                            obj=Material_Interpolation_ISO_SIMP_P3(constParams);
                        otherwise
                            disp('Method not added')
                    end
            end
        end
        computeMatProp()
    end
end