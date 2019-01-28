classdef MemoryManager_MeshUnfitted < MemoryManager
    
    properties (Access = private)
        mesh
        
        subcells
        
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
        
        function saveNewSubcells(obj,subcells,newCellContainingNodes,newCellContainingSubcell)
            obj.subcells = subcells;
            
            obj.updateUpperBounds();
            obj.assignUnfittedNodalProps(newCellContainingNodes);
            obj.assignUnfittedSubcellProps(newCellContainingSubcell);
            obj.assignUnfittedCutCoordIsoPerCell();
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
        
        function updateUpperBounds(obj)
            obj.upperBound_A = obj.lowerBound_A + obj.subcells.nNodes;
            obj.upperBound_B = obj.lowerBound_B + obj.subcells.nSubcells;
            obj.upperBound_C = obj.lowerBound_C + obj.subcells.nSubcells;
        end
        
        function updateLowerBounds(obj)
            obj.lowerBound_A = obj.upperBound_A;
            obj.lowerBound_B = obj.upperBound_B;
            obj.lowerBound_C = obj.upperBound_C;
        end
        
        function assignUnfittedNodalProps(obj,newCellContainingNodes)
            obj.coord_iso(1+obj.lowerBound_A:obj.upperBound_A,:) = obj.subcells.coord_iso;
            obj.coord_global_raw(1+obj.lowerBound_A:obj.upperBound_A,:) = obj.subcells.coord_global;
            obj.cell_containing_nodes(1+obj.lowerBound_A:obj.upperBound_A,:) = newCellContainingNodes;
            obj.levelSet_unfitted(1+obj.lowerBound_A:obj.upperBound_A) = obj.subcells.levelSet;
        end
        
        function assignUnfittedSubcellProps(obj,newCellContainingSubcell)
            obj.connec_local(1+obj.lowerBound_B:obj.upperBound_B,:) = obj.subcells.connec;
            obj.cellContainingSubcell(1+obj.lowerBound_B:obj.upperBound_B,:) = newCellContainingSubcell;
        end
        
        function assignUnfittedCutCoordIsoPerCell(obj)
            for idime = 1:obj.mesh.ndim
                c = obj.subcells.coord_iso(:,idime);
                obj.coord_iso_per_cell(obj.lowerBound_C+1:obj.upperBound_C,:,idime) = c(obj.subcells.connec);
            end
        end
        
    end
    
end

