classdef AllNodesInElemComputer < handle
    
   properties (Access = public)
       allNodesInElem
       cutNodesInElem
   end
   
   properties (Access = private)
   end

   properties (Access = private)
       firstCutEdge
       vertexNodesInElem
       finalNodeNumber
   end
   
   methods (Access = public)
       
       function obj = AllNodesInElemComputer(cParams)
          obj.init(cParams)
       end
       
       function compute(obj)
           obj.computeCutNodePerElem();
           obj.computeAllNodesInElem();
       end
       
   end
    
   methods (Access = private)
       
       function init(obj,cParams)
           obj.firstCutEdge      = cParams.firstCutEdge;
           obj.vertexNodesInElem = cParams.connec;
           obj.finalNodeNumber   = cParams.finalNodeNumber;
       end
       
        function computeAllNodesInElem(obj)
            vertexNodes = obj.vertexNodesInElem;
            cutNodes    = obj.cutNodesInElem;
            allNodes    = [vertexNodes,cutNodes];
            obj.allNodesInElem = allNodes;
        end
        
        function cutNode = computeCutNodePerElem(obj)
            firstCutEdgePerElem = obj.firstCutEdge;
            cutNode = firstCutEdgePerElem + obj.finalNodeNumber;
            obj.cutNodesInElem = cutNode;
        end
       
   end
    
end