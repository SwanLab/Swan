classdef UnfittedMesh_AbstractBuilder < handle
    properties (GetAccess = public, SetAccess = private, Abstract)
        meshType
        maxSubcells
        nnodesSubcell
        
        subcellsMesher
        cutPointsCalculator
        meshPlotter
    end
end

