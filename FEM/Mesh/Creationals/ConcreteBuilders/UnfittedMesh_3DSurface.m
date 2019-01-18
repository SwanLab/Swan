classdef UnfittedMesh_3DSurface < UnfittedMesh_AbstractBuilder
    properties (GetAccess = public, SetAccess = private)
        meshType = 'BOUNDARY';
        max_subcells = 6;
        nnodes_subcell = 3;
        
        subcellsMesher = SubcellsMesher_Boundary_3D;
        cutPointsCalculator = CutPointsCalculator_3D;
        meshPlotter = MeshPlotter_Boundary_3D;
    end
end

