classdef Mesh_Unfitted < Mesh
    properties
        x_fitted
        x_unfitted
        
        full_cells
        empty_cells
        cut_cells
        
        % !! WEIRD THAT IT HAS TO MENTION ITSELF TO AVOID CONFUSION !!
        % (This must be solved in future refactoring. NOT JUST REMOVE THE "unfitted" PREFIX)
        unfitted_coord_iso
        unfitted_coord_iso_per_cell
        unfitted_coord_global
        unfitted_connec_iso
        unfitted_connec_global
        cell_containing_nodes
        cell_containing_subcell
        dvolu_cut
        
        fitted_mesh
        geometryType
    end
    
    properties (Access = protected)
        fitted_geom_interpolation
        
        unfitted_coord_global_raw
        
        max_subcells
        nnodes_subcell
    end
    
    methods
%         function obj = Mesh_Unfitted
%         end
        
%         function obj = create()
%         end
        
        function computeMesh(obj,x_fitted)
            obj.x_fitted = x_fitted;
            obj.findCutCells;
            obj.computeMesh_Delaunay;
        end
        
        function obj = computeMesh_Delaunay(obj)
            [Nodes_n_CutPoints_iso,active_nodes] = obj.findCutPoints_Iso;
            Nodes_n_CutPoints_global = obj.findCutPoints_Global;
            
            obj.allocateMemory_Delaunay;
            
            lowerBound_A = 0; lowerBound_B = 0; lowerBound_C = 0;
            upperBound_A = 0; upperBound_B = 0; upperBound_C = 0;
            
            for icut = 1:length(obj.cut_cells)
                icell = obj.cut_cells(icut);
                cutPoints_iso = Nodes_n_CutPoints_iso(active_nodes(:,:,icut),:,icut);
                cutPoints_global = Nodes_n_CutPoints_global(active_nodes(:,:,icut),:,icut);
                
                [new_unfitted_coord_iso,new_unfitted_coord_global,new_x_unfitted,new_subcell_connec]...
                    = obj.computeSubcells(obj.fitted_mesh.connec(icell,:),cutPoints_iso,cutPoints_global);
                
                number_new_subcells = size(new_subcell_connec,1);
                number_new_coordinates = size(new_unfitted_coord_iso,1);
                
                new_cell_containing_nodes = repmat(icell,[number_new_coordinates 1]);
                
                upperBound_A = lowerBound_A + number_new_coordinates;
                
                obj.unfitted_coord_iso(1+lowerBound_A:upperBound_A,:) = new_unfitted_coord_iso;
                obj.unfitted_coord_global_raw(1+lowerBound_A:upperBound_A,:) = new_unfitted_coord_global;
                obj.cell_containing_nodes(1+lowerBound_A:upperBound_A,:) = new_cell_containing_nodes;
                obj.x_unfitted(1+lowerBound_A:upperBound_A) = new_x_unfitted;
                
                lowerBound_A = upperBound_A;
                
                new_cell_containing_subcell = repmat(icell,[number_new_subcells 1]);
                upperBound_B = lowerBound_B + number_new_subcells;
                obj.unfitted_connec_iso(1+lowerBound_B:upperBound_B,:) = new_subcell_connec;
                obj.cell_containing_subcell(1+lowerBound_B:upperBound_B,:) = new_cell_containing_subcell;
                lowerBound_B = upperBound_B;
                
                upperBound_C = lowerBound_C + number_new_subcells;
                obj.assignUnfittedCutCoordIsoPerCell(new_unfitted_coord_iso,new_subcell_connec,lowerBound_C,upperBound_C);
                lowerBound_C = upperBound_C;
            end
            obj.cleanExtraAllocatedMemory_Delaunay(upperBound_A,upperBound_B,upperBound_C);
        end
        
        function computeGlobalConnectivities(obj)
            obj.unfitted_coord_global =  unique(obj.unfitted_coord_global_raw,'rows','stable');
            obj.unfitted_connec_global = obj.computeFromLocalToGlobalConnectivities(obj.unfitted_connec_global,obj.unfitted_coord_global_raw,obj.unfitted_coord_global,obj.unfitted_connec_iso,obj.cell_containing_nodes,obj.cell_containing_subcell);
        end
        
        function global_connec = computeFromLocalToGlobalConnectivities(obj,global_connec,local_coord,global_coord,local_connec,cell_containing_nodes,cell_containing_subcell)
            for i = 1:size(local_connec,1) % !! VECTORIZE THIS LOOP !!
                icell = cell_containing_subcell(i);
                indexes_in_global_matrix = obj.findCoordinatesIndexesInGlobalCoordinatesMatrix(local_coord(cell_containing_nodes == icell,:),global_coord);
                global_connec(i,:) = indexes_in_global_matrix(local_connec(i,:));
            end
        end
        
        function allocateMemory_Delaunay(obj)
            number_cut_cells = length(obj.cut_cells);
            obj.unfitted_coord_iso = zeros(number_cut_cells*obj.max_subcells*obj.nnodes_subcell,obj.fitted_mesh.ndim);
            obj.unfitted_coord_global_raw = zeros(number_cut_cells*obj.max_subcells*obj.nnodes_subcell,obj.fitted_mesh.ndim);
            obj.x_unfitted = zeros(number_cut_cells*obj.max_subcells*obj.nnodes_subcell,1);
            obj.unfitted_connec_iso = zeros(number_cut_cells*obj.max_subcells,obj.nnodes_subcell);
            obj.unfitted_coord_iso_per_cell = zeros(number_cut_cells*obj.max_subcells,obj.nnodes_subcell,obj.fitted_mesh.ndim);
            obj.unfitted_connec_global = zeros(number_cut_cells*obj.max_subcells,obj.nnodes_subcell);
            obj.cell_containing_subcell = zeros(number_cut_cells*obj.max_subcells*obj.nnodes_subcell,1);
        end
        
        function cleanExtraAllocatedMemory_Delaunay(obj,upperBound_A,upperBound_B,upperBound_C)
            if length(obj.unfitted_coord_iso) > upperBound_A
                obj.unfitted_coord_iso(upperBound_A+1:end,:) = [];
                obj.unfitted_coord_global_raw(upperBound_A+1:end,:) = [];
                obj.x_unfitted(upperBound_A+1:end) = [];
            end
            if length(obj.unfitted_connec_iso) > upperBound_B
                obj.unfitted_connec_iso(upperBound_B+1:end,:) = [];
                obj.unfitted_connec_global(upperBound_B+1:end,:) = [];
                obj.cell_containing_subcell(upperBound_B+1:end) = [];
            end
            
            if length(obj.unfitted_coord_iso_per_cell) > upperBound_C
                obj.unfitted_coord_iso_per_cell(upperBound_C+1:end,:,:) = [];
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
    
    methods %(Access = private)
        function storeFittedMesh(obj,fitted_mesh,fitted_geom_interpolation)
            obj.fitted_mesh = fitted_mesh;
            obj.fitted_geom_interpolation = fitted_geom_interpolation;
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

