classdef Mesh_Unfitted_Properties < handle
    
    properties (GetAccess = public, SetAccess = protected)
        backgroundFullCells
        backgroundEmptyCells
        backgroundCutCells
        
        backgroundGeomInterpolation
        pdim
    end
    
    properties (Access = public, Dependent)
        nCutCells
    end
    
    properties (GetAccess = public, SetAccess = ?UnfittedMesh_AbstractBuilder)
        unfittedType
    end
    
    properties (GetAccess = public, SetAccess = ?MemoryManager_MeshUnfitted)
        coord_iso
        connec_local
        coord_iso_per_cell
        cellContainingSubcell
    end
    
    properties (GetAccess = protected, SetAccess = ?MemoryManager_MeshUnfitted)
        coord_global_raw
        cellContainingNodes
    end
    
    properties (GetAccess = ?MemoryManager_MeshUnfitted, SetAccess = ?UnfittedMesh_AbstractBuilder)
        maxSubcells
        nnodesSubcell
    end
    
    properties (GetAccess = protected, SetAccess = ?UnfittedMesh_AbstractBuilder)
        subcellsMesher
        cutPointsCalculator
        meshPlotter
        cellsClassifier
        memoryManager
    end
    
    methods
        
        function nCutCells = get.nCutCells(obj)
            nCutCells = length(obj.backgroundCutCells);
        end
        
    end
    
end

