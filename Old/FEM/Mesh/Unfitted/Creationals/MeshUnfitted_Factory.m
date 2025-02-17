classdef MeshUnfitted_Factory < handle
    
    methods (Static, Access = public)
        
        function obj = create(type,mesh_background,interpolation_background)
            switch type
                case 'INTERIOR'
                    switch mesh_background.ndim
                        case 2
                            obj = Mesh_Unfitted_2D_Interior(mesh_background.clone,interpolation_background);
                        case 3
                            obj = Mesh_Unfitted_3D_Interior(mesh_background.clone,interpolation_background);
                    end
                case 'BOUNDARY'
                    switch mesh_background.ndim
                        case 2
                            obj = Mesh_Unfitted_2D_Boundary(mesh_background.clone,interpolation_background);
                        case 3
                            obj = Mesh_Unfitted_3D_Boundary(mesh_background.clone,interpolation_background);
                    end
            end
        end
        
    end
    
end

