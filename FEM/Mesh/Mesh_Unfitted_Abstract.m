classdef Mesh_Unfitted_Abstract < handle
    
    properties (GetAccess = public, SetAccess = protected)
        x_background
        x_unfitted
        meshBackground
    end
    
    methods (Access = public, Abstract)
        computeMesh(obj,levelSet)
    end
    
end