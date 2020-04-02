classdef AllEdges2CutEdgesComputer < handle
    
   properties (Access = private) 
        isEdgeCutInElem    
        nElem
        nCutEdgeByElem
   end
    
   methods (Access = public)
       
       function obj = AllEdges2CutEdgesComputer(cParams)
          obj.init(cParams) 
       end
       
       function edge = compute(obj,allEdges)
            allEdges = transpose(allEdges);
            cutEdges = allEdges(obj.isEdgeCutInElem);
            cutEdges = reshape(cutEdges,obj.nCutEdgeByElem,obj.nElem);
            edge = transpose(cutEdges);           
       end
   end
    
   
   methods (Access = private)
       
       function init(obj,cParams)
           obj.isEdgeCutInElem = cParams.isEdgeCutInElem;
           obj.nElem           = cParams.nElem;
           obj.nCutEdgeByElem  = cParams.nCutEdgeByElem;
       end
       
   end
    
    
end