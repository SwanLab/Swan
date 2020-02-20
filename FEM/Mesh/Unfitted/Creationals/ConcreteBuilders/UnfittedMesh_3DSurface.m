classdef UnfittedMesh_3DSurface < UnfittedMesh_AbstractBuilder
    properties (GetAccess = public, SetAccess = private)
        unfittedType = 'BOUNDARY';
        maxSubcells = 6;
        nnodesSubcell = 3;
        
        subcellsMesher      = SubcellsMesher_Boundary_3D;
        cutPointsCalculator = CutPointsCalculator;
        meshPlotter         = MeshPlotter_Boundary_3D;
    end
end

