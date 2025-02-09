classdef CutEdgesComputer < handle
    
    properties (GetAccess = public, SetAccess = private)
        xCutEdgePoint
        isEdgeCut
    end
    
    properties (Access = private)
        isoCoords
    end
    
    properties (Access = private)
        nodesInEdges
        levelSet
    end
    
    methods (Access = public)
        
        function obj = CutEdgesComputer(cParams)
            obj.init(cParams);
        end
        
        function compute(obj)
            obj.computeCutEdges();
            obj.computeXcutEdgePoint();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.levelSet      = cParams.levelSet;
            obj.nodesInEdges  = cParams.nodesInEdges;
            obj.isoCoords     = [-1 1];
        end
        
        function computeCutEdges(obj)
            nodes = obj.nodesInEdges;
            nodes1 = nodes(:,1);
            nodes2 = nodes(:,2);
            ls1 = obj.levelSet(nodes1);
            ls2 = obj.levelSet(nodes2);
            obj.isEdgeCut = xor(ls1<0,ls2<0);
        end
        
        function computeXcutEdgePoint(obj)
            nodes           = obj.nodesInEdges;
            nodesInCutEdges = nodes(obj.isEdgeCut,:);
            node1 = nodesInCutEdges(:,1);
            node2 = nodesInCutEdges(:,2);
            ls1 = obj.levelSet(node1);
            ls2 = obj.levelSet(node2);
            x1  = obj.isoCoords(1);
            x2  = obj.isoCoords(2);
            xCut = x1+ls1.*(x2-x1)./(ls1-ls2);
           obj.xCutEdgePoint = xCut;
        end
        
        
    end
    
end