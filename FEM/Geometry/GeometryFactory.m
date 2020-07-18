classdef GeometryFactory < handle
    
    methods (Access = public, Static)
        
        function obj = create(cParams)
            
            switch cParams.mesh.ndim                
                case 1
                    obj = Geometry_Volumetric(cParams);                                        
                case 2
                    switch cParams.mesh.kFace                        
                        case -1
                            obj = Geometry_Line(cParams);                            
                        case 0
                            obj = Geometry_Volumetric(cParams);
                    end                    
                case 3
                    switch cParams.mesh.kFace                     
                        case -2
                            obj = Geometry_Line(cParams);                                                        
                        case -1
                            obj = Geometry_Surface(cParams);                                                        
                        case 0
                            obj = Geometry_Volumetric(cParams);
                    end                    
            end
            
        end
        
        
    end
    
end