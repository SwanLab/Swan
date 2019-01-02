classdef UnfittedMesh_AbstractBuilder < handle
    properties (GetAccess = public, SetAccess = private, Abstract)
        meshType
        max_subcells
        nnodes_subcell
        
        subcells_Mesher
        cutPoints_Calculator
    end
end

