classdef Mesh_Unfitted_2D_Interior < Mesh_Unfitted_2D & Mesh_Unfitted_Interior
    
    properties
        
    end
    
    methods
        function obj = Mesh_Unfitted_2D_Interior(fitted_mesh,fitted_geom_interpolation)
            obj.storeFittedMesh(fitted_mesh,fitted_geom_interpolation);
        end
    end
end

