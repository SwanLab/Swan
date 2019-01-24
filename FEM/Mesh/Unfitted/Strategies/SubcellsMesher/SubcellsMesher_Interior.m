classdef SubcellsMesher_Interior < SubcellsMesher_Abstract 
    methods (Access = public)
        function [interior_subcell_coord_iso,interior_subcell_coord_global,interior_subcell_x_value,interior_subcell_connec] = computeSubcells(obj,mesh_background,background_geom_interpolation,x_background,background_cell_connec,cutPoints_iso,cutPoints_global)
            interior_subcell_coord_iso = [background_geom_interpolation.pos_nodes; cutPoints_iso];
            interior_subcell_coord_global = [mesh_background.coord(background_cell_connec,:); cutPoints_global];
            interior_subcell_x_value = [x_background(background_cell_connec); zeros(size(cutPoints_iso,1),1)]';
            
            interior_subcell_connec = obj.computeInteriorSubcellsConnectivities(interior_subcell_coord_iso,interior_subcell_x_value);
        end
    end
    
    methods (Access = protected)
        function subcells_connec = computeInteriorSubcellsConnectivities(obj,subcell_coord_iso,subcell_x_value)
            subcells_connec = obj.computeDelaunay(subcell_coord_iso);
            is_interior = all(subcell_x_value(subcells_connec) <= 0,2);
            
            subcells_connec = subcells_connec(is_interior,:);
        end
    end
end

