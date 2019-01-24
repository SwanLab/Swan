classdef Mesh_Unfitted < Mesh ...
                       & Mesh_Unfitted_Abstract...
                       & Mesh_Unfitted_Properties
    
    methods (Access = public)
        
        function obj = Mesh_Unfitted(meshType,meshBackground,interpolation_background)
            obj.build(meshType,meshBackground.ndim);
            obj.init(meshBackground,interpolation_background);
        end
        
        function computeMesh(obj,levelSet_background)
            obj.updateLevelSet(levelSet_background);
            obj.classifyCells();
            if obj.isLevelSetCuttingMesh()
                obj.computeUnfittedMesh();
                obj.computeGlobalConnectivities();
            else
                obj.setCoordinatesAsBackground();
            end
            obj.computeGeometryType();
        end
        
        function mass = computeMass(obj)
            integrator = Integrator.create(obj);
            M2 = integrator.integrateUnfittedMesh(ones(size(obj.levelSet_background)),obj);
            mass = sum(M2);
        end
        
        function plot(obj)
            h = figure;
            obj.add2plot(axes(h));
            light
            axis equal off
            hold off
        end
        
        function add2plot(obj,ax,removedDim,removedDimCoord)
            meshUnfitted = obj.clone();
            if nargin == 4
                meshUnfitted = obj.meshPlotter.patchRemovedDimension(meshUnfitted,removedDim,removedDimCoord);
            end
            obj.meshPlotter.plot(meshUnfitted,ax);
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,meshBackground,backgroundGeomInterpolation)
            obj.ndim = meshBackground.ndim;
            obj.meshBackground = meshBackground;
            obj.backgroundGeomInterpolation = backgroundGeomInterpolation;
        end
        
        function build(obj,meshType,ndim)
            builder = UnfittedMesh_Builder_Factory.create(meshType,ndim);
            builder.build(obj);
        end
        
        function updateLevelSet(obj,levelSet_background)
            obj.levelSet_background = levelSet_background;
        end
        
        function itIs = isLevelSetCuttingMesh(obj)
            itIs = ~isempty(obj.backgroundCutCells);
        end
        
        function computeGlobalConnectivities(obj)
            obj.coord = unique(obj.coord_global_raw,'rows','stable');
            obj.computeFromLocalToGlobalConnectivities;
        end
        
        function classifyCells(obj)
            [F,E,C] = obj.cellsClassifier.classifyCells(obj.levelSet_background,obj.meshBackground.connec);
            obj.backgroundFullCells = F;
            obj.backgroundEmptyCells = E;
            obj.backgroundCutCells = C;
        end
        
        function setCoordinatesAsBackground(obj)
            obj.coord = obj.meshBackground.coord;
        end
        
        function obj = computeUnfittedMesh(obj)
            obj.memoryManager.link(obj);
            obj.memoryManager.allocateMemory();
            
            [cutPoints_iso,cutPoints_global,real_cutPoints] = obj.computeCutPoints();
            
            for icut = 1:obj.nCutCells
                icell = obj.backgroundCutCells(icut);
                
                [new_coord_iso,new_coord_global,new_x_unfitted,new_subcells_connec] = ...
                    obj.computeSubcells(icut,icell,cutPoints_iso,cutPoints_global,real_cutPoints);
                
                nNewSubcells = size(new_subcells_connec,1);
                nNewCoords = size(new_coord_iso,1);
                
                newCellContainingNodes = repmat(icell,[nNewCoords 1]);
                newCellContainingSubcell = repmat(icell,[nNewSubcells 1]);
                
                obj.memoryManager.saveNewSubcells(new_coord_iso,new_coord_global,new_x_unfitted,new_subcells_connec,...
                    newCellContainingNodes,newCellContainingSubcell,nNewSubcells,nNewCoords);
            end
            obj.memoryManager.freeSpareMemory();
            obj.memoryManager.transferData();
        end
        
        function computeFromLocalToGlobalConnectivities(obj)
            for i = 1:size(obj.connec_local,1)
                icell = obj.cellContainingSubcell(i);
                indexes_in_global_matrix = obj.findIndexesOfCoordinatesAinCoordinateMatrixB(obj.coord_global_raw(obj.cell_containing_nodes == icell,:),obj.coord);
                obj.connec(i,:) = indexes_in_global_matrix(obj.connec_local(i,:));
            end
        end
        
        function [new_coord_iso,new_coord_global,new_x_unfitted,new_subcell_connec] = computeSubcells(obj,icut,icell,cutPoints_iso,cutPoints_global,real_cutPoints)
            currentCell_cutPoints_iso = obj.getCurrentCutPoints(cutPoints_iso,real_cutPoints,icut);
            currentCell_cutPoints_global = obj.getCurrentCutPoints(cutPoints_global,real_cutPoints,icut);
            
            [new_coord_iso,new_coord_global,new_x_unfitted,new_subcell_connec]...
                = obj.subcellsMesher.computeSubcells(obj.meshBackground,obj.backgroundGeomInterpolation,obj.levelSet_background,obj.meshBackground.connec(icell,:),currentCell_cutPoints_iso,currentCell_cutPoints_global);
        end
        
        function [cutPoints_iso,cutPoints_global,real_cutPoints] = computeCutPoints(obj)
            obj.cutPointsCalculator.init(obj);
            [cutPoints_iso,real_cutPoints] = obj.cutPointsCalculator.computeCutPoints_Iso();
            cutPoints_global = obj.cutPointsCalculator.computeCutPoints_Global();
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

