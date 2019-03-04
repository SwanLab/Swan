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
        cellsClassifier
        memoryManager
    end
    
    methods (Access = public)
        
        function obj = UnfittedMesh_AbstractBuilder()
            obj.cellsClassifier = CellsClassifier;
            obj.memoryManager = MemoryManager_MeshUnfitted;
        end
        
    end
    
    methods (Access = ?Mesh_Unfitted)
        
        function build(obj,mesh)
            mesh.unfittedType = obj.unfittedType;
            mesh.maxSubcells = obj.maxSubcells;
            mesh.nnodesSubcell = obj.nnodesSubcell;
            mesh.subcellsMesher = obj.subcellsMesher;
            mesh.cutPointsCalculator = obj.cutPointsCalculator;
            mesh.meshPlotter = obj.meshPlotter;
            mesh.cellsClassifier = obj.cellsClassifier;
            mesh.memoryManager = obj.memoryManager;
        end
        
    end
    
end

