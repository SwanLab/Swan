classdef Mesh_Unfitted < Mesh
    properties
        full_cells
        empty_cells
        cut_cells
        
        coord_iso
        connec_local
        coord_iso_per_cell
        cell_containing_subcell
        
        x_fitted
        x_unfitted
        
        fitted_mesh
        geometryType
    end
    
    properties (Access = protected)
        coord_global_raw
        
        cell_containing_nodes
        
        max_subcells
        nnodes_subcell
        
        fitted_geom_interpolation
    end
    
    methods (Access = public)
        function storeFittedMesh(obj,fitted_mesh,fitted_geom_interpolation)
            obj.fitted_mesh = fitted_mesh;
            obj.fitted_geom_interpolation = fitted_geom_interpolation;
        end
        
        function computeMesh(obj,x_fitted)
            obj.x_fitted = x_fitted;
            obj.findCutCells;
            obj.computeMesh_Delaunay;
        end
        
        function computeGlobalConnectivities(obj)
            obj.coord =  unique(obj.coord_global_raw,'rows','stable');
            obj.computeFromLocalToGlobalConnectivities;
        end
    end
    
    methods (Access = private)
        function findCutCells(obj)
            phi_nodes = obj.x_fitted(obj.fitted_mesh.connec);
            phi_case = sum((sign(phi_nodes)<0),2);
            
            obj.full_cells = phi_case == size(obj.fitted_mesh.connec,2);
            obj.empty_cells = phi_case == 0;
            indexes = (1:size(obj.fitted_mesh.connec,1))';
            obj.cut_cells = indexes(~(obj.full_cells | obj.empty_cells));
        end
        
        function obj = computeMesh_Delaunay(obj)
            [Nodes_n_CutPoints_iso,real_cutPoints] = obj.findCutPoints_Iso;
            Nodes_n_CutPoints_global = obj.findCutPoints_Global;
            
            obj.allocateMemory_Delaunay;
            
            lowerBound_A = 0; lowerBound_B = 0; lowerBound_C = 0;
            for icut = 1:length(obj.cut_cells)
                icell = obj.cut_cells(icut);
                currentCell_cutPoints_iso = obj.getCurrentCutPoints(Nodes_n_CutPoints_iso,real_cutPoints,icut);
                currentCell_cutPoints_global = obj.getCurrentCutPoints(Nodes_n_CutPoints_global,real_cutPoints,icut);
                
                [new_coord_iso,new_coord_global,new_x_unfitted,new_subcell_connec]...
                    = obj.computeSubcells(obj.fitted_mesh.connec(icell,:),currentCell_cutPoints_iso,currentCell_cutPoints_global);
                
                number_new_subcells = size(new_subcell_connec,1);
                number_new_coordinates = size(new_coord_iso,1);
                
                new_cell_containing_nodes = repmat(icell,[number_new_coordinates 1]);
                new_cell_containing_subcell = repmat(icell,[number_new_subcells 1]);
                
                [lowerBound_A,lowerBound_B,lowerBound_C] = obj.saveNewSubcells... % !! add / save / store / assign !! (?)
                    (new_coord_iso,new_coord_global,new_x_unfitted,new_subcell_connec,...
                    new_cell_containing_nodes,new_cell_containing_subcell,number_new_subcells,number_new_coordinates,...
                    lowerBound_A,lowerBound_B,lowerBound_C);
            end
            obj.cleanExtraAllocatedMemory_Delaunay(lowerBound_A,lowerBound_B,lowerBound_C);
        end
        
        function computeFromLocalToGlobalConnectivities(obj)
            indexes_in_global_matrix = obj.findIndexesOfCoordinatesAinCoordinateMatrixB(obj.coord_global_raw,obj.coord);
            connec_global_raw = obj.connec_local + repmat(colon(0,2,2*(size(obj.connec_local,1)-1))',[1 size(obj.connec_local,2)]);
            obj.connec = indexes_in_global_matrix(connec_global_raw);
        end
        
        function [lowerBound_A,lowerBound_B,lowerBound_C] = saveNewSubcells...
                (obj,new_coord_iso,new_coord_global,new_x_unfitted,new_subcell_connec,...
                new_cell_containing_nodes,new_cell_containing_subcell,number_new_subcells,number_new_coordinates,...
                lowerBound_A,lowerBound_B,lowerBound_C)
            
            upperBound_A = lowerBound_A + number_new_coordinates;
            obj.assignUnfittedNodalProps(lowerBound_A,upperBound_A,new_coord_iso,new_coord_global,new_x_unfitted,new_cell_containing_nodes);
            lowerBound_A = upperBound_A;
            
            upperBound_B = lowerBound_B + number_new_subcells;
            obj.assignUnfittedSubcellProps(lowerBound_B,upperBound_B,new_subcell_connec,new_cell_containing_subcell);
            lowerBound_B = upperBound_B;
            
            upperBound_C = lowerBound_C + number_new_subcells;
            obj.assignUnfittedCutCoordIsoPerCell(new_coord_iso,new_subcell_connec,lowerBound_C,upperBound_C);
            lowerBound_C = upperBound_C;
        end
        
        function allocateMemory_Delaunay(obj)
            number_cut_cells = length(obj.cut_cells);
            obj.coord_iso = zeros(number_cut_cells*obj.max_subcells*obj.nnodes_subcell,obj.fitted_mesh.ndim);
            obj.coord_global_raw = zeros(number_cut_cells*obj.max_subcells*obj.nnodes_subcell,obj.fitted_mesh.ndim);
            obj.coord_iso_per_cell = zeros(number_cut_cells*obj.max_subcells,obj.nnodes_subcell,obj.fitted_mesh.ndim);
            obj.connec_local = zeros(number_cut_cells*obj.max_subcells,obj.nnodes_subcell);
            obj.connec = zeros(number_cut_cells*obj.max_subcells,obj.nnodes_subcell);
            obj.x_unfitted = zeros(number_cut_cells*obj.max_subcells*obj.nnodes_subcell,1);
            obj.cell_containing_nodes = zeros(number_cut_cells*obj.max_subcells*obj.nnodes_subcell,1);
            obj.cell_containing_subcell = zeros(number_cut_cells*obj.max_subcells*obj.nnodes_subcell,1);
        end
        
        function cleanExtraAllocatedMemory_Delaunay(obj,upperBound_A,upperBound_B,upperBound_C)
            if length(obj.coord_iso) > upperBound_A
                obj.coord_iso(upperBound_A+1:end,:) = [];
                obj.coord_global_raw(upperBound_A+1:end,:) = [];
                obj.cell_containing_nodes(upperBound_A+1:end) = [];
                obj.x_unfitted(upperBound_A+1:end) = [];
            end
            if length(obj.connec_local) > upperBound_B
                obj.connec_local(upperBound_B+1:end,:) = [];
                obj.connec(upperBound_B+1:end,:) = [];
                obj.cell_containing_subcell(upperBound_B+1:end) = [];
            end
            if length(obj.coord_iso_per_cell) > upperBound_C
                obj.coord_iso_per_cell(upperBound_C+1:end,:,:) = [];
            end
        end
        
        function assignUnfittedNodalProps(obj,lowerBound_A,upperBound_A,new_coord_iso,new_coord_global,new_x_unfitted,new_cell_containing_nodes)
            obj.coord_iso(1+lowerBound_A:upperBound_A,:) = new_coord_iso;
            obj.coord_global_raw(1+lowerBound_A:upperBound_A,:) = new_coord_global;
            obj.cell_containing_nodes(1+lowerBound_A:upperBound_A,:) = new_cell_containing_nodes;
            obj.x_unfitted(1+lowerBound_A:upperBound_A) = new_x_unfitted;
        end
        
        function assignUnfittedSubcellProps(obj,lowerBound_B,upperBound_B,new_subcell_connec,new_cell_containing_subcell)
            obj.connec_local(1+lowerBound_B:upperBound_B,:) = new_subcell_connec;
            obj.cell_containing_subcell(1+lowerBound_B:upperBound_B,:) = new_cell_containing_subcell;
        end
    end
    
    methods (Access = private, Static)
        function cutPoints = getCurrentCutPoints(Nodes_n_CutPoints,real_cutPoints,icut)
            cutPoints = Nodes_n_CutPoints(real_cutPoints(:,:,icut),:,icut);
        end
        
        function indexes = findIndexesOfCoordinatesAinCoordinateMatrixB(A,B)
            indexes = zeros(1,size(A,1));
            for inode = 1:size(A,1)
                match = true(size(B,1),1);
                for idime = 1:size(A,2)
                    match = match & B(:,idime) == A(inode,idime);
                end
                indexes(inode) = find(match,1);
            end
        end
    end
end

