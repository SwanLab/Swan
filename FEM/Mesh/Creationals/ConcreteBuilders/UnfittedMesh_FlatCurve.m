classdef UnfittedMesh_FlatCurve < UnfittedMesh_AbstractBuilder
    properties (GetAccess = public, SetAccess = private)
        meshType = 'BOUNDARY';
        max_subcells = 2;
        nnodes_subcell = 2;
        
        subcellsMesher = SubcellsMesher_Boundary_2D;
        cutPointsCalculator = CutPointsCalculator_2D;
        meshPlotter = MeshPlotter_Boundary_2D;
    end
end

