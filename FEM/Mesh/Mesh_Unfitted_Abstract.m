classdef Mesh_Unfitted_Abstract < handle
    
    properties (GetAccess = public, SetAccess = protected)
        levelSet_background
        levelSet_unfitted
        meshBackground
    end
    
    methods (Access = public, Abstract)
        computeMesh(obj,levelSet)
    end
    
end