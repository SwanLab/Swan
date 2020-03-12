classdef MaterialInterpolationFactory < handle
    
    methods (Access = public, Static)
        
        function obj = create(cParams)
            switch cParams.typeOfMaterial
                case 'ISOTROPIC'
                    switch cParams.interpolation
                        case 'SIMPALL'
                            if ~isfield(cParams,'simpAllType') 
                                if isempty(cParams.simpAllType) 
                                    cParams.simpAllType = 'EXPLICIT';
                                end
                            end
                            switch cParams.dim
                                case '2D'
                                    switch cParams.simpAllType
                                        case 'EXPLICIT'
                                            obj = SimpallInterpolationExplicit2D(cParams);
                                        case 'IMPLICIT'
                                            obj = Material_Interpolation_ISO_SIMPALL_2D(cParams);
                                    end
                                case '3D'
                                    switch cParams.simpAllType
                                        case 'EXPLICIT'                                    
                                            obj = Material_Interpolation_ISO_SIMPALL_3D(cParams);
                                        case 'IMPLICIT'
                                            obj = SimpallInterpolationExplicit3D(cParams);
                                    end
                                    
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