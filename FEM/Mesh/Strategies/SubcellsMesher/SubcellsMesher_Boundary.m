classdef SubcellsMesher_Boundary < SubcellsMesher_Abstract
    methods (Access = public) %(Access = ?Mesh_Unfitted)
        function [facets_coord_iso,facets_coord_global,facets_x_value,facets_connec] = computeSubcells(obj,~,background_geom_interpolation,x_background,background_cell_connec,cutPoints_iso,cutPoints_global)
            facets_coord_iso = cutPoints_iso;
            facets_coord_global = cutPoints_global;
            facets_x_value = zeros(1,size(cutPoints_iso,1));
            cell_x_value = x_background(background_cell_connec)';
            
            interior_subcell_coord_iso = [background_geom_interpolation.pos_nodes; cutPoints_iso];
            
            number_nodes = size(background_geom_interpolation.pos_nodes,1);
            facets_connec = obj.computeFacetsConnectivities(facets_coord_iso,interior_subcell_coord_iso,cell_x_value,number_nodes);
        end
    end
    
    methods (Access = public, Abstract) %(Access = ?Mesh_Unfitted, Abstract)
        facets_connec = computeFacetsConnectivities(obj,facets_coord_iso,interior_subcell_coord_iso,cell_x_value)
    end
end

