classdef Mesh_Unfitted_Abstract < handle      
    methods (Access = public, Abstract)
        computeMesh(obj,levelSet)
    end
end