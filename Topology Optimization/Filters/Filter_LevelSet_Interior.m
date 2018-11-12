classdef Filter_LevelSet_Interior < Filter_LevelSet
    properties (Access = public)
        domainType = 'INTERIOR';
    end
    
    %                 function M2 = computeRHS(obj,F1)
    %                     obj.unfitted_mesh.computeDvoluCut;
    %
    %                     posgp_iso = obj.computePosGP(obj.unfitted_mesh.coord_iso_per_cell,obj.interpolation_unfitted,obj.quadrature_unfitted);
    %                     obj.interpolation.computeShapeDeriv(squeeze(posgp_iso));
    %
    %                     shapeValues_CutCells = obj.integrateCutCells(obj.unfitted_mesh.cell_containing_subcell,obj.unfitted_mesh.dvolu_cut);
    %                     shapeValues_All = obj.assembleShapeValues(shapeValues_CutCells,obj.shapeValues_FullCells);
    %
    %                     M2 = obj.rearrangeOutputRHS(shapeValues_All);
    %                 end
    
end