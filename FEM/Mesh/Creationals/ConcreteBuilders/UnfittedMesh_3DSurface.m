classdef UnfittedMesh_3DSurface < UnfittedMesh_AbstractBuilder
    properties (GetAccess = public, SetAccess = private)
        type = 'BOUNDARY';
        max_subcells = 6;
        nnodes_subcell = 3;
        
        subcells_Mesher = SubcellsMesher_Boundary_3D;
        cutPoints_Calculator = CutPoints_Calculator_3D;
    end
end

