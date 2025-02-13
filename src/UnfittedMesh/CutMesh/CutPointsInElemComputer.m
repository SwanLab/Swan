classdef CutPointsInElemComputer < handle
    
    properties (GetAccess = public, SetAccess = private)        
        allNodesInElem
        cutNodesInElem
        xAllNodesInElem
        edgeCutPointInElem
        xCutInElem
        nodesInCutEdges
    end
    
    properties (Access = private)

    end
    
    properties (Access = private)
        cutEdgeInElem
        isEdgeCutInElem
        
        isEdgeCut
        edgesInElem
        nEdgeByElem
        nCutEdges
        nElem

        all2Cut
        
        allNodesinElemParams
        allNodesInElemCoordParams
    end
    
    methods (Access = public)
        
        function obj = CutPointsInElemComputer(cParams)
            obj.init(cParams)
        end
        
        function compute(obj)
            obj.computeEdgeCutPointsInElem();
            obj.computeAllNodesInElem();
            obj.computeXallNodesInElemAndXcut();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.isEdgeCutInElem           = cParams.isEdgeCutInElem;
            obj.isEdgeCut                 = cParams.isEdgeCut;
            obj.edgesInElem               = cParams.edgesInElem;
            obj.nEdgeByElem               = cParams.nEdgeByElem;
            obj.allNodesinElemParams      = cParams.allNodesinElemParams;
            obj.allNodesInElemCoordParams = cParams.allNodesInElemCoordParams;
            obj.all2Cut                   = cParams.all2Cut;
            obj.nElem             = size(obj.edgesInElem,1);
            obj.nCutEdges         = sum(obj.isEdgeCut);
        end
        

        function computeEdgeCutPointsInElem(obj)
            edges = obj.edgesInElem;
            cEdgeInElem = obj.all2Cut.compute(edges);
            nAllEdges = size(obj.isEdgeCut,1);
            cutPoint = zeros(nAllEdges,1);
            cutPoint(obj.isEdgeCut) = 1:obj.nCutEdges;
            edgeCutPoint = zeros(obj.nElem,obj.all2Cut.nCutEdgeByElem);
            for iedge = 1:obj.all2Cut.nCutEdgeByElem
                edge = cEdgeInElem(:,iedge);
                edgeCutPoint(:,iedge) = cutPoint(edge);
            end
            obj.edgeCutPointInElem = edgeCutPoint;
        end
        
        function computeAllNodesInElem(obj)
            edge = obj.edgeCutPointInElem;
            s = obj.allNodesinElemParams;
            s.firstCutEdge = edge;
            aComputer = AllNodesInElemComputer(s);
            aComputer.compute();
            obj.allNodesInElem = aComputer.allNodesInElem;
            obj.cutNodesInElem = aComputer.cutNodesInElem;
        end
        
        function computeXallNodesInElemAndXcut(obj)
            s = obj.allNodesInElemCoordParams;
            s.all2Cut               = obj.all2Cut;
            s.edgeCutPointInElem    = obj.edgeCutPointInElem;
            a = AllNodesInElemCoordinatesComputer(s);
            a.compute();
            obj.xAllNodesInElem = a.xAllNodesInElem;
            obj.xCutInElem = a.xCutInElem;
            obj.nodesInCutEdges = a.nodesInCutEdges;
        end
        
    end
    
end