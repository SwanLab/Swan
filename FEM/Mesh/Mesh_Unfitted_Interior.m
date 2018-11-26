classdef Mesh_Unfitted_Interior < Mesh_Unfitted % !! Change to Interior Contour / Skin !!
    properties %(GetAccess = public, SetAccess = protected)
        dvolu_cut
    end
    
    methods (Access = public)
        function [interior_subcell_coord_iso,interior_subcell_coord_global,interior_subcell_x_value,interior_subcell_connec] = computeSubcells(obj,background_cell_connec,cutPoints_iso,cutPoints_global)
            interior_subcell_coord_iso = [obj.background_geom_interpolation.pos_nodes; cutPoints_iso];
            interior_subcell_coord_global = [obj.mesh_background.coord(background_cell_connec,:); cutPoints_global];
            interior_subcell_x_value = [obj.x_background(background_cell_connec); zeros(size(cutPoints_iso,1),1)]';
            
            interior_subcell_connec = obj.computeInteriorSubcellsConnectivities(interior_subcell_coord_iso,interior_subcell_x_value);
        end
    end
    
    methods (Static, Access = private)
        function subcell_connec = computeInteriorSubcellsConnectivities(subcell_coord_iso,subcell_x_value)
            DT = delaunayTriangulation(subcell_coord_iso);
            subcell_connec = DT.ConnectivityList;
            is_interior = all(subcell_x_value(subcell_connec) <= 0,2);
            
            subcell_connec = subcell_connec(is_interior,:);
        end
    end
end

