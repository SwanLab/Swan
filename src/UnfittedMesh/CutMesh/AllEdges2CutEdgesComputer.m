classdef AllEdges2CutEdgesComputer < handle
    
    properties (GetAccess = public, SetAccess = private)
        nCutEdgeByElem
    end
    
   properties (Access = private)
        isEdgeCutInElem
   end
    
   methods (Access = public)
       
       function obj = AllEdges2CutEdgesComputer(cParams)
          obj.init(cParams)
       end
       
       function edge = compute(obj,allEdges)
            nElem = size(allEdges,1);
            allEdges = transpose(allEdges);
            cutEdges = allEdges(obj.isEdgeCutInElem);
            cutEdges = reshape(cutEdges,obj.nCutEdgeByElem,nElem);
            edge     = transpose(cutEdges);
       end
       
       
   end
    
   
   methods (Access = private)
       
       function init(obj,cParams)
           obj.isEdgeCutInElem = cParams.isEdgeCutInElem;
           obj.nCutEdgeByElem = unique(sum(obj.isEdgeCutInElem,1));
       end
       
   end
    
    
end