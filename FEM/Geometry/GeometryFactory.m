classdef GeometryFactory < handle
    
    methods (Access = public, Static)
        
        function obj = create(cParams)
            switch cParams.mesh.geometryType
                case 'Line'
                    obj = Geometry_Line(cParams);
                case 'Surface'
                    switch cParams.mesh.ndim
                        case 2
                            obj = Geometry_Volumetric(cParams);
                        case 3
                            obj = Geometry_Surface(cParams);
                    end
                case 'Volume'
                    obj = Geometry_Volumetric(cParams);
            end
        end
    end
    
end