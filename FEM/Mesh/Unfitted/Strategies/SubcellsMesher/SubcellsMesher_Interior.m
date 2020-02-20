classdef SubcellsMesher_Interior < SubcellsMesher_Abstract
    
    properties (Access = protected)
        interiorSubcells
        subcells_connec
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
        
        function computeAllPossibleSubcellsInCell(obj)
            if size(obj.coord_iso,2) == 1
                [~,I] = sort(obj.coord_iso);
                connec = [I circshift(I,-1)];
                connec(end,:) = [];
                obj.subcells_connec = connec;
            else
                obj.subcells_connec = obj.computeDelaunay(obj.coord_iso);
            end
        end
        
    end
    
    methods (Access = private)
        
        function computeCoordsIso(obj)
            nodes = obj.mesh.backgroundGeomInterpolation.pos_nodes;
            cutPts = obj.cutPoints.ISO;
            obj.coord_iso = obj.patch(nodes,cutPts);
        end
        
        function computeCoordsGlobal(obj)
            nodes = obj.mesh.meshBackground.coord(obj.cell_connec,:);
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

