classdef UnfittedMesh_FlatCurve < UnfittedMesh_AbstractBuilder
    properties (GetAccess = public, SetAccess = private)
        unfittedType = 'BOUNDARY';
        maxSubcells = 2;
        nnodesSubcell = 2;
        
        subcellsMesher      = SubcellsMesher_Boundary_2D;
        cutPointsCalculator = CutPointsCalculator;
        meshPlotter         = MeshPlotter_Boundary_2D;
    end
end

