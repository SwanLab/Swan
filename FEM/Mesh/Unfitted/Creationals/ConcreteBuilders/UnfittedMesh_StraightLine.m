classdef UnfittedMesh_StraightLine < UnfittedMesh_AbstractBuilder
    properties (GetAccess = public, SetAccess = private)
        unfittedType = 'INTERIOR';
        maxSubcells = 2;
        nnodesSubcell = 2;
        
        subcellsMesher      = SubcellsMesher_1D;
        cutPointsCalculator = CutPointsCalculator;
        meshPlotter         = MeshPlotter_1D;
    end
end

