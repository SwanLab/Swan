classdef Mesh_Unfitted_Interior < Mesh_Unfitted
    
    properties
        
    end
    
    methods
        function obj = Mesh_Unfitted_Interior
        end
        
        function [subcell_coord_iso,subcell_coord_global,subcell_x_value,interior_subcell_connec] = computeSubcells(obj,fitted_cell_connec,subcell_cutPoints_iso,subcell_cutPoints_global)
            subcell_coord_iso = [obj.fitted_geom_interpolation.pos_nodes; subcell_cutPoints_iso];
            subcell_coord_global = [obj.fitted_mesh.coord(fitted_cell_connec,:); subcell_cutPoints_global];
            subcell_x_value = [obj.x_fitted(fitted_cell_connec); zeros(size(subcell_cutPoints_iso,1),1)]';
            
            interior_subcell_connec = obj.computeInteriorSubcellsConnectivities(subcell_coord_iso,subcell_x_value);
        end
    end
end

