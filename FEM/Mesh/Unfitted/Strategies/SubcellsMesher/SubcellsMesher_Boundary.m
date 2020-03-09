classdef SubcellsMesher_Boundary < SubcellsMesher
    
    properties (GetAccess = protected, SetAccess = private)
        nCellNodes
    end
    
    methods (Access = protected, Abstract)
        
        computeFacetsConnectivities(obj)
        
    end
    
    methods (Access = protected)
        
        function computeCoordinates(obj)
            obj.coord_iso = obj.cutPoints.ISO;
            obj.coord_global = obj.cutPoints.GLOBAL;
        end
        
        function computeConnectivities(obj)
            obj.computeFacetsConnectivities();
        end
        
        function computeLevelSet(obj)
            obj.levelSet = obj.boundary_levelSet;
        end
        
    end
    
    methods
        
        function nCellNodes = get.nCellNodes(obj)
            nCellNodes = size(obj.posNodes,1);
        end
        
    end
end

