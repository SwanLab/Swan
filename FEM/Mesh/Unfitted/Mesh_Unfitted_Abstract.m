classdef Mesh_Unfitted_Abstract < handle
    
    properties (GetAccess = public, SetAccess = protected)
        meshBackground
        levelSet_background
        levelSet_unfitted
    end
    
    methods (Access = public, Abstract)
        
        computeMesh(obj,levelSet)
        
    end
    
    methods (Access = public)
        
        function setLevelSetUnfitted(obj,LS)
            obj.levelSet_unfitted = LS;
        end
        
        function setLevelSetBackground(obj,LS)
            obj.levelSet_background = LS;
        end
        
    end
    
end