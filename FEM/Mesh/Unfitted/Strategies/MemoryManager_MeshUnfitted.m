classdef MemoryManager_MeshUnfitted < MemoryManager
    
    properties (Access = private)
        mesh
        
        connec
        levelSet_unfitted
        coord_iso
        connec_local
        coord_iso_per_cell
        cellContainingSubcell
        coord_global_raw
        cell_containing_nodes
        
        lowerBound_A
        lowerBound_B
        lowerBound_C
        upperBound_A
        upperBound_B
        upperBound_C
    end
    
    methods (Access = public)
        
        function obj = MemoryManager_MeshUnfitted()
            obj.init();
        end
        
        function link(obj,mesh)
            obj.mesh = mesh;
        end
        
        function allocateMemory(obj)
            obj.coord_iso = zeros(obj.mesh.nCutCells*obj.mesh.maxSubcells*obj.mesh.nnodesSubcell,obj.mesh.meshBackground.ndim);
            obj.coord_global_raw = zeros(obj.mesh.nCutCells*obj.mesh.maxSubcells*obj.mesh.nnodesSubcell,obj.mesh.meshBackground.ndim);
            obj.coord_iso_per_cell = zeros(obj.mesh.nCutCells*obj.mesh.maxSubcells,obj.mesh.nnodesSubcell,obj.mesh.meshBackground.ndim);
            obj.connec_local = zeros(obj.mesh.nCutCells*obj.mesh.maxSubcells,obj.mesh.nnodesSubcell);
            obj.connec = (zeros(obj.mesh.nCutCells*obj.mesh.maxSubcells,obj.mesh.nnodesSubcell));
            obj.levelSet_unfitted = (zeros(obj.mesh.nCutCells*obj.mesh.maxSubcells*obj.mesh.nnodesSubcell,1));
            obj.cell_containing_nodes = zeros(obj.mesh.nCutCells*obj.mesh.maxSubcells*obj.mesh.nnodesSubcell,1);
            obj.cellContainingSubcell = zeros(obj.mesh.nCutCells*obj.mesh.maxSubcells*obj.mesh.nnodesSubcell,1);
        end
        
        function freeSpareMemory(obj)
            if length(obj.coord_iso) > obj.upperBound_A
                obj.coord_iso(obj.upperBound_A+1:end,:) = [];
                obj.coord_global_raw(obj.upperBound_A+1:end,:) = [];
                obj.cell_containing_nodes(obj.upperBound_A+1:end) = [];
                obj.levelSet_unfitted(obj.upperBound_A+1:end) = [];
            end
            
            if length(obj.connec_local) > obj.upperBound_B
                obj.connec_local(obj.upperBound_B+1:end,:) = [];
                obj.connec(obj.upperBound_B+1:end,:) = [];
                obj.cellContainingSubcell(obj.upperBound_B+1:end) = [];
            end
            
            if length(obj.coord_iso_per_cell) > obj.upperBound_C
                obj.coord_iso_per_cell(obj.upperBound_C+1:end,:,:) = [];
            end
        end
        
        function saveNewSubcells(obj,new_coord_iso,new_coord_global,newLevelSet_unfitted,new_subcell_connec,...
                newCellContainingNodes,newCellContainingSubcell,nNewSubcells,nNewCoords)
            
            obj.updateUpperBounds(nNewCoords,nNewSubcells);
            obj.assignUnfittedNodalProps(new_coord_iso,new_coord_global,newLevelSet_unfitted,newCellContainingNodes);
            obj.assignUnfittedSubcellProps(new_subcell_connec,newCellContainingSubcell);
            obj.assignUnfittedCutCoordIsoPerCell(new_coord_iso,new_subcell_connec);
            obj.updateLowerBounds();
        end
        
        function transferData(obj)
            obj.mesh.coord_iso = obj.coord_iso;
            obj.mesh.coord_global_raw = obj.coord_global_raw;
            obj.mesh.coord_iso_per_cell = obj.coord_iso_per_cell;
            obj.mesh.connec_local = obj.connec_local;
            obj.mesh.setConnec(obj.connec);
            obj.mesh.setLevelSetUnfitted(obj.levelSet_unfitted);
            obj.mesh.cell_containing_nodes = obj.cell_containing_nodes;
            obj.mesh.cellContainingSubcell = obj.cellContainingSubcell;
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.lowerBound_A = 0;
            obj.lowerBound_B = 0;
            obj.lowerBound_C = 0;
        end
        
        function updateUpperBounds(obj,nNewCoords,nNewSubcells)
            obj.upperBound_A = obj.lowerBound_A + nNewCoords;
            obj.upperBound_B = obj.lowerBound_B + nNewSubcells;
            obj.upperBound_C = obj.lowerBound_C + nNewSubcells;
        end
        
        function updateLowerBounds(obj)
            obj.lowerBound_A = obj.upperBound_A;
            obj.lowerBound_B = obj.upperBound_B;
            obj.lowerBound_C = obj.upperBound_C;
        end
        
        function assignUnfittedNodalProps(obj,new_coord_iso,new_coord_global,new_x_unfitted,newCellContainingNodes)
            obj.coord_iso(1+obj.lowerBound_A:obj.upperBound_A,:) = new_coord_iso;
            obj.coord_global_raw(1+obj.lowerBound_A:obj.upperBound_A,:) = new_coord_global;
            obj.cell_containing_nodes(1+obj.lowerBound_A:obj.upperBound_A,:) = newCellContainingNodes;
            obj.levelSet_unfitted(1+obj.lowerBound_A:obj.upperBound_A) = new_x_unfitted;
        end
        
        function assignUnfittedSubcellProps(obj,new_subcell_connec,newCellContainingSubcell)
            obj.connec_local(1+obj.lowerBound_B:obj.upperBound_B,:) = new_subcell_connec;
            obj.cellContainingSubcell(1+obj.lowerBound_B:obj.upperBound_B,:) = newCellContainingSubcell;
        end
        
        function assignUnfittedCutCoordIsoPerCell(obj,new_coord_iso,new_interior_subcell_connec)
            for idime = 1:obj.mesh.ndim
                new_coord_iso_ = new_coord_iso(:,idime);
                obj.coord_iso_per_cell(obj.lowerBound_C+1:obj.upperBound_C,:,idime) = new_coord_iso_(new_interior_subcell_connec);
            end
        end
        
    end
    
end

