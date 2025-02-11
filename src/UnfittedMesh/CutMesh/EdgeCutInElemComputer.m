classdef EdgeCutInElemComputer < handle
    
    properties (Access = private)
       nElem
       nEdgeByElem
    end
    
    properties (Access = private)
      edgesInElem
      isEdgeCut
    end
    
    methods (Access = public)
        
        function obj = EdgeCutInElemComputer(cParams)
            obj.init(cParams)
        end
        
        function isEdgeCut = compute(obj)
            isEdgeCut = false(obj.nEdgeByElem,obj.nElem);
            for iedge = 1:obj.nEdgeByElem
                edge = obj.edgesInElem(:,iedge);
                isEdgeCut(iedge,:) = obj.isEdgeCut(edge);
            end
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.edgesInElem = cParams.edgesInElem;
            obj.isEdgeCut   = cParams.isEdgeCut;
            obj.nElem       = size(obj.edgesInElem,1);
            obj.nEdgeByElem = size(obj.edgesInElem,2);
        end
        
    end
    
end