classdef GeometryFactory < handle
    
    methods (Access = public, Static)
        
        function obj = create(cParams)
            nGeom = cParams.mesh.ndim + cParams.mesh.kFace; 
            switch nGeom
                case 1
                   obj = Geometry_Line(cParams);                            
                case 2
                    switch cParams.mesh.ndim
                        case 2
                            obj = Geometry_Volumetric(cParams);
                        case 3
                            obj = Geometry_Surface(cParams);                                                        
                    end
                case 3
                    obj = Geometry_Volumetric(cParams);
            end
        end
     end
    
end