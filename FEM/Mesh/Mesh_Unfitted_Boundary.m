classdef Mesh_Unfitted_Boundary < Mesh_Unfitted
    methods
        function [facets_coord_iso,facets_coord_global,facets_x_value,facets_connec] = computeSubcells(obj,fitted_cell_connec,cutPoints_iso,cutPoints_global)
            facets_coord_iso = cutPoints_iso;
            facets_coord_global = cutPoints_global;
            facets_x_value = zeros(1,size(cutPoints_iso,1));
            cell_x_value = obj.x_fitted(fitted_cell_connec)';
            
            interior_subcell_coord_iso = [obj.fitted_geom_interpolation.pos_nodes; cutPoints_iso];
            
            facets_connec = obj.computeFacetsConnectivities(facets_coord_iso,interior_subcell_coord_iso,cell_x_value);
        end
    end
end

