classdef Mesh_Unfitted_2D_Interior < Mesh_Unfitted_2D & Mesh_Unfitted_Interior
    methods
        function obj = Mesh_Unfitted_2D_Interior(fitted_mesh,fitted_geom_interpolation)
            obj.storeFittedMesh(fitted_mesh,fitted_geom_interpolation);
            obj.geometryType = 'TRIANGLE';
            obj.max_subcells = 6;
            obj.nnodes_subcell = 3;
        end
        
        function computeDvoluCut(obj)
            x1 = obj.coord_iso_per_cell(:,1,1);   y1 = obj.coord_iso_per_cell(:,1,2);
            x2 = obj.coord_iso_per_cell(:,2,1);   y2 = obj.coord_iso_per_cell(:,2,2);
            x3 = obj.coord_iso_per_cell(:,3,1);   y3 = obj.coord_iso_per_cell(:,3,2);
            obj.dvolu_cut = 0.5*abs((x2-x1).*(y3-y1)-(x3-x1).*(y2-y1));
        end
    end
end

