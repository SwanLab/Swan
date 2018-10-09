classdef Mesh_Unfitted_3D_Boundary < Mesh_Unfitted_3D & Mesh_Unfitted_Boundary
    methods
        function obj = Mesh_Unfitted_3D_Boundary(fitted_mesh,fitted_geom_interpolation)
            obj.storeFittedMesh(fitted_mesh,fitted_geom_interpolation);
            obj.geometryType = 'TRIANGLE';
            obj.max_subcells = 6; % !! ?? !!
            obj.nnodes_subcell = 3;
        end
    end
end
