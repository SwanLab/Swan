classdef Mesh_Unfitted_3D_Interior < Mesh_Unfitted_3D & Mesh_Unfitted_Interior
    methods
        function obj = Mesh_Unfitted_3D_Interior(fitted_mesh,fitted_geom_interpolation)
            obj.storeFittedMesh(fitted_mesh,fitted_geom_interpolation);
            obj.geometryType = 'TETRAHEDRA';
            obj.max_subcells = 20;
            obj.nnodes_subcell = 4;
        end
    end
end

