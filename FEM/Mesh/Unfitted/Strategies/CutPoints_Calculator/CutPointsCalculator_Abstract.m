classdef CutPointsCalculator_Abstract < handle
    
    properties(GetAccess = public, SetAccess = protected)
        cutPointsIso
        cutPointsGlobal
        activeCutPoints
    end
    
    properties (GetAccess = protected, SetAccess = private)
        meshBackground
        levelSet_background
        backgroundCutCells
        backgroundGeomInterpolation
    end
    
    methods (Access = protected, Abstract)
        
        computeCutPoints_Iso(obj)
        computeCutPoints_Global(obj)
        
    end
    
    methods (Access = public)
        
        function init(obj,mesh)
            obj.meshBackground = mesh.meshBackground;
            obj.levelSet_background = mesh.levelSet_background;
            obj.backgroundCutCells = mesh.backgroundCutCells;
            obj.backgroundGeomInterpolation = mesh.backgroundGeomInterpolation;
        end
        
        function computeCutPoints(obj)
            obj.computeCutPoints_Iso();
            obj.computeCutPoints_Global();
        end
        
        function cutPoints = getThisCellCutPoints(obj,i)
            cutPoints.ISO = obj.getThisCellActiveCutPointsIso(i);
            cutPoints.GLOBAL = obj.getThisCellActiveCutPointsGlobal(i);
        end
        
    end
    
    methods (Access = private)
        
        function cP = getThisCellActiveCutPointsIso(obj,i)
            cP = obj.cutPointsIso(obj.activeCutPoints(:,:,i),:,i);
        end
        
        function cP = getThisCellActiveCutPointsGlobal(obj,i)
            cP = obj.cutPointsGlobal(obj.activeCutPoints(:,:,i),:,i);
        end
        
        function cP = getThisCellActiveCutPoints(obj,cutPoints,i)
            cP = cutPoints(obj.activeCutPoints(:,:,i),:,i);
        end
        
    end
    
end

