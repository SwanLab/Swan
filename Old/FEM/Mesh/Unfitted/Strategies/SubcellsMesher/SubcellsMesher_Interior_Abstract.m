classdef SubcellsMesher_Interior_Abstract < SubcellsMesher
    
    properties (Access = protected)
        interiorSubcells
        subcells_connec
    end
    
    methods (Access = protected, Abstract)
        computeAllPossibleSubcellsInCell(obj)
    end
    
    methods (Access = protected)
        
        function computeCoordinates(obj)
            obj.computeCoordsIso();
            obj.computeCoordsGlobal();
        end
        
        function computeLevelSet(obj)
            obj.levelSet = obj.patch(obj.cell_levelSet,obj.boundary_levelSet)';
        end
        
        function computeConnectivities(obj)
            obj.computeInteriorSubcellsConnectivities();
        end
        
        function computeInteriorSubcellsConnectivities(obj)
            obj.computeAllPossibleSubcellsInCell();
            obj.getInteriorSubcells();
        end
        
    end
    
    methods (Access = private)
        
        function computeCoordsIso(obj)
            nodes = obj.posNodes;
            cutPts = obj.cutPoints.ISO;
            obj.coord_iso = obj.patch(nodes,cutPts);
        end
        
        function computeCoordsGlobal(obj)
            nodes = obj.coordsBackground(obj.cell_connec,:);
            cutPts = obj.cutPoints.GLOBAL;
            obj.coord_global = obj.patch(nodes,cutPts);
        end
        
        function getInteriorSubcells(obj)
            obj.findInteriorSubcells();
            obj.connec = obj.subcells_connec(obj.interiorSubcells,:);
        end
        
        function findInteriorSubcells(obj)
            isNodeInterior = obj.levelSet(obj.subcells_connec) <= 0;
            obj.interiorSubcells = all(isNodeInterior,2);
        end
        
    end
    
end

