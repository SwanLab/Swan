classdef UnfittedMesh_Volumetric < UnfittedMesh_AbstractBuilder
    properties (GetAccess = public, SetAccess = private)
        type = 'INTERIOR';
        max_subcells = 20;
        nnodes_subcell = 4;
        
        subcells_Mesher = SubcellsMesher_Interior;
        cutPoints_Calculator = CutPoints_Calculator_3D;
    end
end

