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
    end
    
    properties (Access = protected)
        fitted_geom_interpolation
        
        unfitted_cut_coord_global_raw
        
        max_subcells
        nnodes_subcell
    end
    
    methods
        function obj = Mesh_Unfitted(fitted_mesh,x_fitted,fitted_geom_interpolation)
            obj.fitted_mesh = fitted_mesh;
            obj.fitted_geom_interpolation = fitted_geom_interpolation;
            obj.x_fitted = x_fitted;
        end
        
        function computeCutMesh(obj)
            obj.findCutCells;
            obj.computeCutMesh_Delaunay;
        end
        
        function obj = computeCutMesh_Delaunay(obj)
            [nodes_n_cutpoints_iso,active_nodes] = obj.findCutPoints_Iso;
            nodes_n_cutpoints_global = obj.findCutPoints_Global;
            
            obj.computeDelaunay_allocateMemory;
            
            k0 = 0; m0 = 0; c0 = 0;
            k1 = 0; m1 = 0; c1 = 0;
            for icut = 1:length(obj.cut_cells)
                icell = obj.cut_cells(icut);
                subcell_cutPoints_iso = nodes_n_cutpoints_iso(active_nodes(:,:,icut),:,icut);
                subcell_cutPoints_global = nodes_n_cutpoints_global(active_nodes(:,:,icut),:,icut);
                
                [new_unfitted_cut_coord_iso,new_unfitted_cut_coord_global,new_x_unfitted_cut,new_subcell_cut_interior_connec_iso]...
                    = obj.computeInteriorSubcells(obj.fitted_mesh.connec(icell,:),subcell_cutPoints_iso,subcell_cutPoints_global);
                
                new_nodes_containing_cell = repmat(icell,[size(new_unfitted_cut_coord_iso,1) 1]);
                
                k1 = k0 + size(new_unfitted_cut_coord_iso,1);
                obj.unfitted_cut_coord_iso(1+k0:k1,:) = new_unfitted_cut_coord_iso;
                obj.unfitted_cut_coord_global_raw(1+k0:k1,:) = new_unfitted_cut_coord_global;
                obj.nodes_containing_cell(1+k0:k1,:) = new_nodes_containing_cell;
                obj.x_unfitted_cut(1+k0:k1) = new_x_unfitted_cut;
                k0 = k1;
                
                new_subcell_containing_cell = repmat(icell,[size(new_subcell_cut_interior_connec_iso,1) 1]);
                m1 = m0+size(new_subcell_cut_interior_connec_iso,1);
                obj.unfitted_cut_connec_iso(1+m0:m1,:) = new_subcell_cut_interior_connec_iso;
                obj.subcell_containing_cell(1+m0:m1,:) = new_subcell_containing_cell;
                m0 = m1;
                
                c1 = c0 + size(new_subcell_cut_interior_connec_iso,1);
                obj.assignUnfittedCutCoordIsoPerCell(new_unfitted_cut_coord_iso,new_subcell_cut_interior_connec_iso,c0,c1);
                c0 = c1;
            end
            
            obj.cleanExtraAllocatedMemory(k1,m1,c1);
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
        
        function computeDelaunay_allocateMemory(obj)
            obj.unfitted_cut_coord_iso = zeros(length(obj.cut_cells)*obj.max_subcells*obj.nnodes_subcell,obj.fitted_mesh.ndim);
            obj.unfitted_cut_coord_global_raw = zeros(length(obj.cut_cells)*obj.max_subcells*obj.nnodes_subcell,obj.fitted_mesh.ndim);
            obj.x_unfitted_cut = zeros(length(obj.cut_cells)*obj.max_subcells*obj.nnodes_subcell,1);
            obj.unfitted_cut_connec_iso = zeros(length(obj.cut_cells)*obj.max_subcells,obj.nnodes_subcell);
            obj.unfitted_cut_coord_iso_per_cell = zeros(length(obj.cut_cells)*obj.max_subcells,obj.nnodes_subcell,obj.fitted_mesh.ndim);
            obj.unfitted_cut_connec_global = zeros(length(obj.cut_cells)*obj.max_subcells,obj.nnodes_subcell);
            obj.subcell_containing_cell = zeros(length(obj.cut_cells)*obj.max_subcells*obj.nnodes_subcell,1);
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
        
        function subcell_connec = computeInteriorSubcellsConnectivities(subcell_coord_iso,subcell_x_value)
            DT = delaunayTriangulation(subcell_coord_iso);
            subcell_connec = DT.ConnectivityList;
            is_interior = all(subcell_x_value(subcell_connec) <= 0,2);
            
            subcell_connec = subcell_connec(is_interior,:);
        end
    end
end

