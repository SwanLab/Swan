classdef SubcellsMesher < handle
    
    properties (GetAccess = public, SetAccess = protected)
        subcells
    end
    
    properties (GetAccess = protected, SetAccess = protected)
        coord_iso
        coord_global
        levelSet
        connec
    end
    
    properties (Access = private)
    end
    
    properties (GetAccess = protected, SetAccess = private)
        mesh
        interior_coord_iso
        cutPoints
        nCutPoints
        cell_nodes
        cell_connec
        cell_levelSet
        boundary_levelSet
        coordsBackground
        posNodes
        levelSetBackground

    end
    
    methods (Access = protected, Abstract)
        
        computeCoordinates(obj)
        computeConnectivities(obj)
        computeLevelSet(obj)
        
    end
    
    methods (Access = public, Static)
        
        function obj = create(cParams)
            f = SubcellsMesherFactory();
            obj = f.create(cParams);
        end
        
    end
    
    methods (Access = public)
 
        function computeSubcells(obj,cParams)

            cellConnec = cParams.cellConnec;
            cutPoints  = cParams.cutPoints;
      
            obj.saveInputs(cellConnec,cutPoints);
            obj.computeLevelSet();
            obj.computeCoordinates();
            obj.computeConnectivities();
            obj.packResultsInStruct();
        end
        
    end
    
    methods (Access = protected)
        
        function init(obj,cParams)
            obj.posNodes = cParams.posNodes;
            obj.levelSetBackground = cParams.levelSetBackground;
            obj.coordsBackground = cParams.coordsBackground;
        end 
    end
    
    methods (Access = protected, Static)
        
        function connectivities = computeDelaunay(coordinates)
            DT = delaunayTriangulation(coordinates);
            connectivities = DT.ConnectivityList;
        end
        
        function w = patch(u,v)
            w = [u;v];
        end
        
    end
    
    methods (Access = private)
        
        function saveInputs(obj,cell_connec,cutPoints)
            obj.cell_connec = cell_connec;
            obj.cutPoints = cutPoints;
        end
        
        function packResultsInStruct(obj)
            obj.subcells.levelSet     = obj.levelSet;
            obj.subcells.coord_iso    = obj.coord_iso;
            obj.subcells.coord_global = obj.coord_global;
            obj.subcells.connec       = obj.connec;
            
            obj.subcells.nNodes       = size(obj.coord_iso,1);
            obj.subcells.nSubcells    = size(obj.connec,1);
        end
        
    end
    
    methods
        
        function coord = get.interior_coord_iso(obj)
            coord = obj.patch(obj.cell_nodes,obj.cutPoints.ISO);
        end
        
        function n = get.nCutPoints(obj)
            n = size(obj.cutPoints.ISO,1);
        end
        
        function nodes = get.cell_nodes(obj)
            nodes = obj.posNodes;
        end
        
        function lvlSet = get.cell_levelSet(obj)
            lvlSet = obj.levelSetBackground(obj.cell_connec);
        end
        
        function lvlSet = get.boundary_levelSet(obj)
            lvlSet = zeros(obj.nCutPoints,1);
        end
        
    end
end

