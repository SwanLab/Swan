classdef MaterialInterpolationFactory < handle
    
    methods (Access = public, Static)
        
        function obj = create(cParams)            
            switch cParams.dim 
                case '2D'
                  cParams.ndim = 2;                                    
                case '3D'
                  cParams.ndim = 3;                                    
            end
               
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
                                            obj = SimpAllInterpolationExplicit2D(cParams);
                                        case 'IMPLICIT'
                                            obj = SimpAllInterpolationImplicit2D(cParams);
                                    end
                                case '3D'
                                    switch cParams.simpAllType
                                        case 'EXPLICIT'                                    
                                            obj = SimpAllInterpolationExplicit3D(cParams);
                                        case 'IMPLICIT'
                                            obj = SimpAllInterpolationImplicit3D(cParams);
                                    end
                                    
                            end
                        case 'SIMP_Adaptative'
                            obj = SimpInterpolationAdaptative(cParams);
                        case 'SIMP_P3'
                            obj = SimpInterpolationP3(cParams);
                        otherwise
                            error('Invalid Material Interpolation method.');
                    end
                otherwise
                    error('Invalid type of material');
            end
            
        end
        
        
    end
    
    
    
end