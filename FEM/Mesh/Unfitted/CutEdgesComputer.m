classdef CutEdgesComputer < handle
    
    properties (GetAccess = public, SetAccess = private)

        nodesInCutEdges
      
        cutEdges        
        elemCases        
        firstCutEdge
    end
    
    properties (Access = private)
        edgesComputer
        levelSet
        isEdgeCutInElem
        nElem
        nCutEdges
        cutEdgeInElem
        nCutEdgeByElem  
        code
    end
    
    methods (Access = public)
        
        function obj = CutEdgesComputer(cParams)
            obj.init(cParams);
        end
        
        function compute(obj)
            obj.computeCutEdges();
            obj.computeNcutEdges();
            obj.computeNodesInCutEdges();
            obj.computeIsEdgeCutInElem();
            obj.computeCutEdgeInElem();
            obj.computeElementCases();
            obj.computeCutNodePerElemen();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.levelSet      = cParams.levelSet;
            obj.edgesComputer = cParams.edgesComputer;
            obj.nElem = size(obj.edgesComputer.edgesInElem,1);
            obj.nCutEdgeByElem = 2;
            obj.code = [5 3 6];
        end
        
        function computeCutEdges(obj)
            nodes = obj.edgesComputer.nodesInEdges;
            nodes1 = nodes(:,1);
            nodes2 = nodes(:,2);
            ls1 = obj.levelSet(nodes1);
            ls2 = obj.levelSet(nodes2);
            obj.cutEdges = xor(ls1<0,ls2<0);
        end
        
        function computeNodesInCutEdges(obj)
            nodes = obj.edgesComputer.nodesInEdges;
            cEdges = obj.cutEdges;
            obj.nodesInCutEdges = nodes(cEdges,:);            
        end
        
        function computeNcutEdges(obj)
            obj.nCutEdges = sum(obj.cutEdges);
        end
        
        function isEdgeCut = computeIsEdgeCutInElem(obj)
            edgesInElem = obj.edgesComputer.edgesInElem;
            nEdgeByElem = obj.edgesComputer.nEdgeByElem;
            isEdgeCut = false(nEdgeByElem,obj.nElem);
            for iedge = 1:nEdgeByElem
                edge = edgesInElem(:,iedge);
                isEdgeCut(iedge,:) = obj.cutEdges(edge);
            end
            obj.isEdgeCutInElem = isEdgeCut;
        end
        
        function computeCutEdgeInElem(obj)
            edgesInCutElemens = obj.edgesComputer.edgesInElem;            
            edgesInCutElemens = transpose(edgesInCutElemens);
            edges = edgesInCutElemens(obj.isEdgeCutInElem);
            edges = reshape(edges,obj.nCutEdgeByElem,obj.nElem);
            edge = transpose(edges);
            obj.cutEdgeInElem = edge;
        end
                
        function computeElementCases(obj)
            edgeCases(:,1) = obj.computeEdgeCases();
            nCases = size(unique(edgeCases),1);
            obj.elemCases = false(obj.nElem,nCases);
            for icase = 1:nCases
                isEdgeCase = edgeCases == obj.code(icase);
                obj.elemCases(:,icase) = isEdgeCase;
            end
        end
        
        function d = computeEdgeCases(obj)
            isCut = obj.isEdgeCutInElem;
            nEdgeByElem = obj.edgesComputer.nEdgeByElem;
            edges = (1:nEdgeByElem) - 1;
            pow2vector = 2.^(edges);
            d = pow2vector*isCut;
        end        
        
        function computeCutNodePerElemen(obj)
            nAllEdges = size(obj.cutEdges,1);
            cEdges = zeros(nAllEdges,1);
            cEdges(obj.cutEdges) = 1:obj.nCutEdges;
            cutEdge = zeros(obj.nElem,obj.nCutEdgeByElem);            
            for iedge = 1:obj.nCutEdgeByElem
                edge = obj.cutEdgeInElem(:,iedge);                
                cutEdge(:,iedge) = cEdges(edge);
            end           
            obj.firstCutEdge = cutEdge;
        end
        
    end
    
end