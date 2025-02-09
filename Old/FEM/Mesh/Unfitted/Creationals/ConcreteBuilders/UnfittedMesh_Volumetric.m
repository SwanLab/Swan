classdef UnfittedMesh_Volumetric < UnfittedMesh_AbstractBuilder
    properties (GetAccess = public, SetAccess = private)
        unfittedType = 'INTERIOR';
        maxSubcells = 20;
        nnodesSubcell = 4;
        
        subcellsMesher      = SubcellsMesher_Interior;
        cutPointsCalculator = CutPointsCalculator;
        meshPlotter         = MeshPlotter_Null;
    end
end

