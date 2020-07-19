classdef MemoryManager_MeshUnfitted < MemoryManager
    
    properties (GetAccess = public, SetAccess = private)
        connec
        levelSet_unfitted
        coord_iso
        connec_local
        subcellIsoCoords
        cellContainingSubcell
        coord_global_raw
        cellContainingNodes        
    end
    
    properties (Access = private)
        mesh
        
        subcells
        
        lowerBound_A
        lowerBound_B
        lowerBound_C
        upperBound_A
        upperBound_B
        upperBound_C
        
        ndimIso
        unfittedType
        
        maxSubcells
        nnodesSubCell
        
        nCutCells
        ndim
    end
    
    methods (Access = public)
        
        function obj = MemoryManager_MeshUnfitted(cParams)
            obj.unfittedType = cParams.unfittedType;
            obj.ndimIso      = cParams.ndimIso;
            obj.nCutCells    = cParams.nCutCells;
            obj.ndim         = cParams.ndim;            
            obj.init();            
        end
              
        function allocateMemory(obj)
            obj.init();
            
            nCell  = obj.nCutCells*obj.maxSubcells;
            nNodes = obj.nCutCells*obj.maxSubcells*obj.nnodesSubCell;
            
            obj.coord_iso             = zeros(nNodes,obj.ndim);
            obj.coord_global_raw      = zeros(nNodes,obj.ndim);
            obj.subcellIsoCoords      = zeros(nCell,obj.nnodesSubCell,obj.ndim);
            obj.connec_local          = zeros(nCell,obj.nnodesSubCell);
            obj.connec                = zeros(nCell,obj.nnodesSubCell);
            obj.levelSet_unfitted     = zeros(nNodes,1);
            obj.cellContainingNodes   = zeros(nNodes,1);
            obj.cellContainingSubcell = zeros(nNodes,1);
        end
        
        function freeSpareMemory(obj)
            if length(obj.coord_iso) > obj.upperBound_A
                obj.coord_iso(obj.upperBound_A+1:end,:) = [];
                obj.coord_global_raw(obj.upperBound_A+1:end,:) = [];
                obj.cellContainingNodes(obj.upperBound_A+1:end) = [];
                obj.levelSet_unfitted(obj.upperBound_A+1:end) = [];
            end
            
            if length(obj.connec_local) > obj.upperBound_B
                obj.connec_local(obj.upperBound_B+1:end,:) = [];
                obj.connec(obj.upperBound_B+1:end,:) = [];
                obj.cellContainingSubcell(obj.upperBound_B+1:end) = [];
            end
            
            if length(obj.subcellIsoCoords) > obj.upperBound_C
                obj.subcellIsoCoords(obj.upperBound_C+1:end,:,:) = [];
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
        
  
         function link(obj,mesh)
%             obj.mesh = mesh;
         end        
        
         function transferData(obj)
         end
        
    end
    
    methods (Access = private)
     
        function init(obj)
            obj.lowerBound_A = 0;
            obj.lowerBound_B = 0;
            obj.lowerBound_C = 0;
            obj.computeMaxAndNodesSubcells();            
        end
        
        function computeMaxAndNodesSubcells(obj)
            switch obj.ndimIso
                case 'Line'
                    obj.maxSubcells = 2;
                    obj.nnodesSubCell = 2;
                case 'Surface'
                    switch obj.unfittedType
                        case 'INTERIOR'
                            obj.maxSubcells = 6;
                            obj.nnodesSubCell = 3;
                        case 'BOUNDARY'
                            obj.maxSubcells = 2;
                            obj.nnodesSubCell = 2;
                    end
                case 'Volume'
                    switch obj.unfittedType
                        case 'INTERIOR'
                            obj.maxSubcells = 20;
                            obj.nnodesSubCell = 4;
                        case 'BOUNDARY'
                            obj.maxSubcells = 6;
                            obj.nnodesSubCell = 3;
                    end
            end
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
            indexA = 1+obj.lowerBound_A:obj.upperBound_A;
            obj.coord_iso(indexA,:) = obj.subcells.coord_iso;
            obj.coord_global_raw(indexA,:) = obj.subcells.coord_global;
            obj.cellContainingNodes(indexA,:) = newCellContainingNodes;
            obj.levelSet_unfitted(indexA) = obj.subcells.levelSet;
        end
        
        function assignUnfittedSubcellProps(obj,newCellContainingSubcell)
            indexB = 1+obj.lowerBound_B:obj.upperBound_B;
            obj.connec_local(indexB,:)          = obj.subcells.connec;
            obj.cellContainingSubcell(indexB,:) = newCellContainingSubcell;
        end
        
        function assignUnfittedCutCoordIsoPerCell(obj)
            for idime = 1:obj.ndim
                c = obj.subcells.coord_iso(:,idime);
                indexC = obj.lowerBound_C+1:obj.upperBound_C;
                obj.subcellIsoCoords(indexC,:,idime) = c(obj.subcells.connec);
            end
        end
        
    end
    
end

