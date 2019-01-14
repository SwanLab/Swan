classdef Mesh_Unfitted < Mesh & Mesh_Unfitted_Abstract
    properties (Access = public) %(Access = private) Assigned by a Builder
        subcells_Mesher
        cutPoints_Calculator
    end
    
    properties (GetAccess = public, SetAccess = private)
        pdim
        meshType
        
        max_subcells
        nnodes_subcell
        
        coord_iso
        connec_local
        coord_iso_per_cell
        cell_containing_subcell
        
        background_full_cells
        background_empty_cells
        background_cut_cells
        
        x_background
        x_unfitted
        
        meshBackground
    end
    
    %    properties (GetAccess = public, SetAccess = protected) % !! Change to private?? !!
    %         integrationParams = struct(...
    %             coord_iso...
    %             connec_local...
    %             coord_iso_per_cell...
    %             cell_containing_subcell...
    %             background_full_cells...
    %             background_empty_cells...
    %             background_cut_cells...
    %             x_background...
    %             x_unfitted...
    %             meshBackground...
    %             );
    %     end
    
    properties (Access = protected)
        coord_global_raw
        cell_containing_nodes
        background_geom_interpolation
    end
    
    methods (Access = public)
        function obj = Mesh_Unfitted(meshType,meshBackground,interpolation_background)
            obj.build(meshType,meshBackground.ndim);
            obj.init(meshBackground,interpolation_background);
        end
        
        function computeMesh(obj,x_background)
            obj.x_background = x_background;
            obj.findCutCells;
            if ~isempty(obj.background_cut_cells)
                obj.computeMesh_Delaunay;
            else
                obj.coord = obj.meshBackground.coord;
                obj.connec = obj.computeDelaunay(obj.coord);
                
            end
            obj.computeGlobalConnectivities;
            obj.computeGeometryType;
        end
        
        function plot(obj)
            h = figure;
            obj.add2plot(axes(h));
            light
            axis equal off
            hold off
        end
        
        function mass = computeMass(obj)
            integrator = Integrator.create(obj);
            M2 = integrator.integrateUnfittedMesh(ones(size(obj.x_background)),obj);
            mass = sum(M2);
        end
    end
    
    methods (Access = private)
        function init(obj,meshBackground,background_geom_interpolation)
            obj.ndim = meshBackground.ndim;
            obj.meshBackground = meshBackground;
            obj.background_geom_interpolation = background_geom_interpolation;
        end
        
        function build(obj,meshType,ndim)
            builder = UnfittedMesh_Builder_Factory.create(meshType,ndim);
            obj.meshType = builder.meshType;
            obj.max_subcells = builder.max_subcells;
            obj.nnodes_subcell = builder.nnodes_subcell;
            obj.subcells_Mesher =	builder.subcells_Mesher;
            obj.cutPoints_Calculator = builder.cutPoints_Calculator;
        end
        
        function computeGlobalConnectivities(obj)
            obj.coord =  unique(obj.coord_global_raw,'rows','stable');
            obj.computeFromLocalToGlobalConnectivities;
        end
        
        function findCutCells(obj)
            phi_nodes = obj.x_background(obj.meshBackground.connec);
            phi_case = sum((sign(phi_nodes)<0),2);
            
            obj.background_full_cells = phi_case == size(obj.meshBackground.connec,2);
            obj.background_empty_cells = phi_case == 0;
            indexes = (1:size(obj.meshBackground.connec,1))';
            obj.background_cut_cells = indexes(~(obj.background_full_cells | obj.background_empty_cells));
        end
        
        function obj = computeMesh_Delaunay(obj)
            [Nodes_n_CutPoints_iso,real_cutPoints] = obj.computeCutPoints_Iso;
            Nodes_n_CutPoints_global = obj.computeCutPoints_Global;
            
            obj.allocateMemory_Delaunay();
            
            lowerBound_A = 0; lowerBound_B = 0; lowerBound_C = 0;
            for icut = 1:length(obj.background_cut_cells)
                icell = obj.background_cut_cells(icut);
                
                [new_coord_iso,new_coord_global,new_x_unfitted,new_subcell_connec] = obj.computeSubcells(icut,icell,Nodes_n_CutPoints_iso,Nodes_n_CutPoints_global,real_cutPoints);
                
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
            %             nnode = size(obj.connec_local,2);
            %             indexes_in_global_matrix = obj.findIndexesOfCoordinatesAinCoordinateMatrixB(obj.coord_global_raw,obj.coord);
            %             connec_global_raw = obj.connec_local + repmat(colon(0,nnode,nnode*(size(obj.connec_local,1)-1))',[1 nnode]);
            %             obj.connec = indexes_in_global_matrix(connec_global_raw);
            
            for i = 1:size(obj.connec_local,1) % !! VECTORIZE THIS LOOP !!
                icell = obj.cell_containing_subcell(i);
                indexes_in_global_matrix = obj.findIndexesOfCoordinatesAinCoordinateMatrixB(obj.coord_global_raw(obj.cell_containing_nodes == icell,:),obj.coord);
                obj.connec(i,:) = indexes_in_global_matrix(obj.connec_local(i,:));
            end
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
            obj.assignUnfittedCutCoordIsoPerCell(lowerBound_C,upperBound_C,new_coord_iso,new_subcell_connec);
            lowerBound_C = upperBound_C;
        end
        
        function allocateMemory_Delaunay(obj)
            number_cut_cells = length(obj.background_cut_cells);
            obj.coord_iso = zeros(number_cut_cells*obj.max_subcells*obj.nnodes_subcell,obj.meshBackground.ndim);
            obj.coord_global_raw = zeros(number_cut_cells*obj.max_subcells*obj.nnodes_subcell,obj.meshBackground.ndim);
            obj.coord_iso_per_cell = zeros(number_cut_cells*obj.max_subcells,obj.nnodes_subcell,obj.meshBackground.ndim);
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
        
        function assignUnfittedCutCoordIsoPerCell(obj,lowerBound_C,upperBound_C,new_coord_iso,new_interior_subcell_connec)
            for idime = 1:obj.ndim
                new_coord_iso_ = new_coord_iso(:,idime);
                obj.coord_iso_per_cell(lowerBound_C+1:upperBound_C,:,idime) = new_coord_iso_(new_interior_subcell_connec);
            end
        end
        
        
        function [new_coord_iso,new_coord_global,new_x_unfitted,new_subcell_connec] = computeSubcells(obj,icut,icell,Nodes_n_CutPoints_iso,Nodes_n_CutPoints_global,real_cutPoints)
            currentCell_cutPoints_iso = obj.getCurrentCutPoints(Nodes_n_CutPoints_iso,real_cutPoints,icut);
            currentCell_cutPoints_global = obj.getCurrentCutPoints(Nodes_n_CutPoints_global,real_cutPoints,icut);
            [new_coord_iso,new_coord_global,new_x_unfitted,new_subcell_connec]...
                = obj.subcells_Mesher.computeSubcells(obj.meshBackground,obj.background_geom_interpolation,obj.x_background,obj.meshBackground.connec(icell,:),currentCell_cutPoints_iso,currentCell_cutPoints_global);
        end
        
        function [Nodes_n_CutPoints_iso,real_cutPoints] = computeCutPoints_Iso(obj)
            [Nodes_n_CutPoints_iso,real_cutPoints] = obj.cutPoints_Calculator.computeCutPoints_Iso(obj.meshBackground,obj.x_background,obj.background_cut_cells,obj.background_geom_interpolation);
        end
        
        function [Nodes_n_CutPoints_global,real_cutPoints] = computeCutPoints_Global(obj)
            [Nodes_n_CutPoints_global,real_cutPoints] = obj.cutPoints_Calculator.computeCutPoints_Global(obj.meshBackground,obj.x_background,obj.background_cut_cells,obj.background_geom_interpolation);
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

