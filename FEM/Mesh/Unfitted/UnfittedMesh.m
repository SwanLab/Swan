classdef UnfittedMesh < handle
    
    properties (GetAccess = public, SetAccess = private)
        unfittedType
        meshBackground
        levelSet_background
        backgroundCutCells
        geometryType
        coord
        connec
        coord_iso_per_cell
        cellContainingSubcell
        backgroundFullCells
        
        nActiveMeshes
    end
    
    properties (Access = private)
        oldUnfittedMesh
    end
    
    
    methods (Access = public)
        
        function obj = UnfittedMesh(cParams)
            obj.oldUnfittedMesh = Mesh_Unfitted.create2(cParams);
        end
        
        function compute(obj,lvlSet)
            obj.oldUnfittedMesh.computeMesh(lvlSet);
        end
        
        function m = computeMass(obj)
            m = obj.oldUnfittedMesh.computeMass();
        end
        
        function aMeshes = getActiveMeshes(obj)
            aMeshes = obj.oldUnfittedMesh.getActiveMeshes();
        end
        
    end
    
    methods
        
        function type = get.unfittedType(obj)
            type = obj.oldUnfittedMesh.unfittedType;
        end
        
        function mB = get.meshBackground(obj)
            mB = obj.oldUnfittedMesh.meshBackground;
        end
        
        function lvlSet = get.levelSet_background(obj)
            lvlSet = obj.oldUnfittedMesh.levelSet_background;
        end
        
        function backCutCells = get.backgroundCutCells(obj)
            backCutCells = obj.oldUnfittedMesh.backgroundCutCells;
        end
        
        function gType = get.geometryType(obj)
            gType = obj.oldUnfittedMesh.geometryType;
        end
        
        function coord = get.coord(obj)
            coord = obj.oldUnfittedMesh.coord;
        end
        
        function connec = get.connec(obj)
            connec = obj.oldUnfittedMesh.connec;
        end
        
        function coord_iso_per_cell = get.coord_iso_per_cell(obj)
            coord_iso_per_cell = obj.oldUnfittedMesh.coord_iso_per_cell;
        end
        
        function cellContainingSubcell = get.cellContainingSubcell(obj)
            cellContainingSubcell = obj.oldUnfittedMesh.cellContainingSubcell;
        end
        
        function backgroundFullCells = get.backgroundFullCells(obj)
            backgroundFullCells = obj.oldUnfittedMesh.backgroundFullCells;
        end
        
        function nActiveMeshes = get.nActiveMeshes(obj)
            nActiveMeshes = obj.oldUnfittedMesh.nActiveMeshes;
        end
        
    end
    
end