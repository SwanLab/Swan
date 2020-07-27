classdef CutPointsInElemComputer < handle
    
    properties (GetAccess = public, SetAccess = private)        
        allNodesInElem
        xAllNodesInElem
        edgeCutPointInElem
        xCut
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
        nCutEdgeByElem

        all2Cut
        
        allNodesinElemParams
        allNodesInElemCoordParams        
    end
    
    methods (Access = public)
        
        function obj = CutPointsInElemComputer(cParams)
            obj.init(cParams)
        end
        
        function compute(obj)
            obj.computeIsEdgeCutInElem();
            obj.computeNcutEdgeByElem();
            obj.createAllEdges2CutEdge();
            obj.computeCutEdgeInElem();
            obj.computeEdgeCutPointsInElem();
            obj.computeAllNodesInElem();
            obj.computeXallNodesInElemAndXcut();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.isEdgeCut                 = cParams.isEdgeCut;
            obj.edgesInElem               = cParams.edgesInElem;
            obj.nEdgeByElem               = cParams.nEdgeByElem;
            obj.allNodesinElemParams      = cParams.allNodesinElemParams;
            obj.allNodesInElemCoordParams = cParams.allNodesInElemCoordParams;
            obj.nElem             = size(obj.edgesInElem,1);
            obj.nCutEdges         = sum(obj.isEdgeCut);
            obj.nCutEdgeByElem    = 2;
        end
        
        function isEdgeCut = computeIsEdgeCutInElem(obj)
            isEdgeCut = false(obj.nEdgeByElem,obj.nElem);
            for iedge = 1:obj.nEdgeByElem
                edge = obj.edgesInElem(:,iedge);
                isEdgeCut(iedge,:) = obj.isEdgeCut(edge);
            end
            obj.isEdgeCutInElem = isEdgeCut;
        end
        
        function computeNcutEdgeByElem(obj)
            n = unique(sum(obj.isEdgeCutInElem));
            obj.nCutEdgeByElem = n;            
        end
        
        function createAllEdges2CutEdge(obj)
           s.isEdgeCutInElem = obj.isEdgeCutInElem;
           s.nElem           = obj.nElem;
           s.nCutEdgeByElem  = obj.nCutEdgeByElem;
           obj.all2Cut = AllEdges2CutEdgesComputer(s);
        end
        
        function computeCutEdgeInElem(obj)   
           edges = obj.edgesInElem;
           obj.cutEdgeInElem = obj.all2Cut.compute(edges);
        end
        
        function computeEdgeCutPointsInElem(obj)
            cEdgeInElem = obj.cutEdgeInElem;
            nAllEdges = size(obj.isEdgeCut,1);
            cutPoint = zeros(nAllEdges,1);
            cutPoint(obj.isEdgeCut) = 1:obj.nCutEdges;
            edgeCutPoint = zeros(obj.nElem,obj.nCutEdgeByElem);
            for iedge = 1:obj.nCutEdgeByElem
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
        end
        
        function computeXallNodesInElemAndXcut(obj)
            s = obj.allNodesInElemCoordParams;
            s.nCutEdgeByElem        = obj.nCutEdgeByElem;
            s.all2Cut               = obj.all2Cut;        
            s.edgeCutPointInElem    = obj.edgeCutPointInElem;
            a = AllNodesInElemCoordinatesComputer(s);
            a.compute();
            obj.xAllNodesInElem = a.xAllNodesInElem;
            obj.xCut = a.xCut;
        end
        
    end
    
end