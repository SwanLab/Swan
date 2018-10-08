classdef Mesh_Unfitted_3D_Interior < Mesh_Unfitted_3D & Mesh_Unfitted_Interior
    
    properties
        
    end
    
    methods
        function obj = Mesh_Unfitted_3D_Interior(fitted_mesh,fitted_geom_interpolation)
            obj.storeFittedMesh(fitted_mesh,fitted_geom_interpolation);
        end
    end
end

