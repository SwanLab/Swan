classdef CutPointsCalculator_Abstract < handle
    
    properties (GetAccess = protected, SetAccess = private)
        meshBackground
        levelSet_background
        backgroundCutCells
        backgroundGeomInterpolation
    end
    
    methods (Access = public, Abstract)
        
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
        
    end
    
end

