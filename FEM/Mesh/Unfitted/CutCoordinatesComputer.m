classdef CutCoordinatesComputer < handle
    
    properties (GetAccess = public, SetAccess = private)
      coord
    end
    
    properties (Access = private)
        levelSet
        backgroundCoord
        xCutPoints
        nodesInCutEdges
    end
    
    methods (Access = public)
        
        function obj = CutCoordinatesComputer(cParams)
            obj.init(cParams)
        end
        
        function compute(obj)
            obj.computeCutPoints();
            obj.computeCutMeshCoordinates();
        end
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.nodesInCutEdges = cParams.nodesInCutEdges;
            obj.levelSet        = cParams.levelSet;
            obj.backgroundCoord = cParams.backgroundCoord;
        end
        
        function computeCutPoints(obj)
            node1 = obj.nodesInCutEdges(:,1);
            node2 = obj.nodesInCutEdges(:,2);
            ls1 = obj.levelSet(node1);
            ls2 = obj.levelSet(node2);
            x1  = obj.backgroundCoord(node1,:);
            x2  = obj.backgroundCoord(node2,:);
            xCut = x1+ls1.*(x2-x1)./(ls1-ls2);
            obj.xCutPoints = xCut;
        end        
        
        function computeCutMeshCoordinates(obj)
            obj.coord = [obj.backgroundCoord;obj.xCutPoints];
        end
        
    end
    
end