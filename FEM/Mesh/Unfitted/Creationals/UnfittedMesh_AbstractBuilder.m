classdef UnfittedMesh_AbstractBuilder < handle
    
    properties (GetAccess = public, SetAccess = private, Abstract)
        unfittedType
        maxSubcells
        nnodesSubcell
        
        subcellsMesher
        cutPointsCalculator
        meshPlotter
    end
    
    properties (GetAccess = public, SetAccess = private)
      %  cellsClassifier
        memoryManager
    end
    
    methods (Access = public)
        
        function obj = UnfittedMesh_AbstractBuilder()
       %     obj.cellsClassifier = CellsClassifier;
            obj.memoryManager = MemoryManager_MeshUnfitted;
        end
        
    end

    
end

