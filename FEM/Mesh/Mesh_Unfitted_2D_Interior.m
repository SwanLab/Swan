classdef Mesh_Unfitted_2D_Interior < Mesh_Unfitted_2D & Mesh_Unfitted_Interior
    methods
        function obj = Mesh_Unfitted_2D_Interior(fitted_mesh,fitted_geom_interpolation)
            obj.storeFittedMesh(fitted_mesh,fitted_geom_interpolation);
            obj.geometryType = 'TRIANGLE';
            obj.max_subcells = 6;
            obj.nnodes_subcell = 3;
        end
    end
end

