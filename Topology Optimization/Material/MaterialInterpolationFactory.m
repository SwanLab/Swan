classdef MaterialInterpolationFactory < handle
    
    
    methods (Access = public, Static)
        
        function obj = create(cParams)
            switch cParams.typeOfMaterial
                case 'ISOTROPIC'
                    switch cParams.interpolation
                        case 'SIMPALL'
                            switch cParams.dim
                                case '2D'
                                    obj = Material_Interpolation_ISO_SIMPALL_2D(cParams);
                                case '3D'
                                    obj = Material_Interpolation_ISO_SIMPALL_3D(cParams);
                            end
                        case 'SIMP_Adaptative'
                            obj = Material_Interpolation_ISO_SIMP_Adaptative(cParams);
                        case 'SIMP_P3'
                            obj = Material_Interpolation_ISO_SIMP_P3(cParams);
                        otherwise
                            error('Invalid Material Interpolation method.');
                    end
                otherwise
                    error('Invalid type of material');
            end
            
        end
        
        
    end
    
    
    
end