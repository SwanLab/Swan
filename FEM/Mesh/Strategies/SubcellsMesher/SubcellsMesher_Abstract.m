classdef SubcellsMesher_Abstract < handle
    methods (Access = public, Abstract)
        [coord_iso,coord_global,x_unfitted,subcell_connec]...
            = computeSubcells(obj,mesh_background,background_geom_interpolation,x_background,mesh_background_connec,currentCell_cutPoints_iso,currentCell_cutPoints_global)
    end
    
    methods (Access = protected, Static)
        function connectivities = computeDelaunay(coordinates)
            DT = delaunayTriangulation(coordinates);
            connectivities = DT.ConnectivityList;
        end
    end
end

