classdef UnfittedMesh_FlatCurve < UnfittedMesh_AbstractBuilder
    properties (GetAccess = public, SetAccess = private)
        meshType = 'BOUNDARY';
        max_subcells = 2;
        nnodes_subcell = 2;
        
        subcells_Mesher = SubcellsMesher_Boundary_2D;
        cutPoints_Calculator = CutPoints_Calculator_2D;
    end
end

