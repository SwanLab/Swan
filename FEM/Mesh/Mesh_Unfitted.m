classdef Mesh_Unfitted < Mesh
    properties
        x_fitted
        x_unfitted
        x_unfitted_cut
        
        full_cells
        empty_cells
        cut_cells
        
        unfitted_cut_coord_iso
        unfitted_cut_coord_iso_per_cell
        unfitted_cut_coord_global
        unfitted_cut_connec_iso
        unfitted_cut_connec_global
        nodes_containing_cell
        subcell_containing_cell
        dvolu_cut
        
        fitted_mesh
        
        i_debug
        flag_debug=false
        debug_elems
    end
    
    properties (Access = protected)
        fitted_geom_interpolation
        
        unfitted_cut_coord_global_raw
        
        max_subcells
        nnodes_subcell
    end
    
    methods
        function obj = Mesh_Unfitted(fitted_mesh,x_fitted,fitted_geom_interpolation)
            if nargin > 0
                obj.fitted_mesh = fitted_mesh;
                obj.fitted_geom_interpolation = fitted_geom_interpolation;
                obj.x_fitted = x_fitted;
            end
        end
        function computeGlobalConnectivities(obj)
            obj.unfitted_cut_coord_global =  unique(obj.unfitted_cut_coord_global_raw,'rows','stable');
            obj.unfitted_cut_connec_global = obj.computeFromLocalToGlobalConnectivities(obj.unfitted_cut_connec_global,obj.unfitted_cut_coord_global_raw,obj.unfitted_cut_coord_global,obj.unfitted_cut_connec_iso,obj.nodes_containing_cell,obj.subcell_containing_cell);
        end
        
        function [subcell_coord_iso,subcell_coord_global,subcell_x_value,subcell_cut_interior_connec] = computeInteriorSubcells(obj,fitted_cell_connec,subcell_cutPoints_iso,subcell_cutPoints_global)
            subcell_coord_iso = [obj.fitted_geom_interpolation.pos_nodes; subcell_cutPoints_iso];
            subcell_coord_global = [obj.fitted_mesh.coord(fitted_cell_connec,:); subcell_cutPoints_global];
            subcell_x_value = [obj.x_fitted(fitted_cell_connec); zeros(size(subcell_cutPoints_iso,1),1)]';
            
            subcell_cut_interior_connec = obj.computeInteriorSubcellsConnectivities(subcell_coord_iso,subcell_x_value);
        end
        
        function global_connectivities = computeFromLocalToGlobalConnectivities(obj,global_connectivities,local_matrix_coord,global_matrix_coord,local_connec,nodes_containing_cell,subcell_containing_cell)
            for i = 1:size(local_connec,1)
                icell = subcell_containing_cell(i);
                indexes_in_global_matrix = obj.findCoordinatesIndexesInGlobalCoordinatesMatrix(local_matrix_coord(nodes_containing_cell == icell,:),global_matrix_coord);
                global_connectivities(i,:) = indexes_in_global_matrix(local_connec(i,:));
            end
        end
        function cleanExtraAllocatedMemory(obj,k,m,c)
            if length(obj.unfitted_cut_coord_iso) > k
                obj.unfitted_cut_coord_iso(k+1:end,:) = [];
                obj.unfitted_cut_coord_global_raw(k+1:end,:) = [];
                obj.x_unfitted_cut(k+1:end) = [];
            end
            if length(obj.unfitted_cut_connec_iso) > m
                obj.unfitted_cut_connec_iso(m+1:end,:) = [];
                obj.unfitted_cut_connec_global(m+1:end,:) = [];
                obj.subcell_containing_cell(m+1:end) = [];
            end
            
            if length(obj.unfitted_cut_coord_iso_per_cell) > c
                obj.unfitted_cut_coord_iso_per_cell(c+1:end,:,:) = [];
            end
        end
        
        function findCutCells(obj)
            phi_nodes = obj.x_fitted(obj.fitted_mesh.connec);
            phi_case = sum((sign(phi_nodes)<0),2);
            
            obj.full_cells = phi_case == size(obj.fitted_mesh.connec,2);
            obj.empty_cells = phi_case == 0;
            indexes = (1:size(obj.fitted_mesh.connec,1))';
            obj.cut_cells = indexes(~(obj.full_cells | obj.empty_cells));
        end

    end
    
    methods (Static)
        function indexes_in_global_matrix = findCoordinatesIndexesInGlobalCoordinatesMatrix(coordinates_local,coordinates_global)
            indexes_in_global_matrix = zeros(1,size(coordinates_local,1));
            for inode = 1:size(coordinates_local,1)
                match = true(size(coordinates_global,1),1);
                for idime = 1:size(coordinates_local,2)
                    match = match & coordinates_global(:,idime) == coordinates_local(inode,idime);
                end
                indexes_in_global_matrix(inode) = find(match,1);
            end
        end
    end
end

