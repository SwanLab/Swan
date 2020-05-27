classdef SubCellsCasesComputer < handle
    
   properties (GetAccess = public, SetAccess = private)
       subCellCases
   end
   
   properties (Access = private)
       isEdgeCutInElem
       code
   end
   
   methods (Access = public)
       
       function obj = SubCellsCasesComputer(cParams)
          obj.init(cParams) 
       end
       
       function compute(obj)
            nElem  = size(obj.isEdgeCutInElem,2);
            nCases = length(obj.code);
            edgeCases(:,1) = obj.computeEdgeCases();            
            obj.subCellCases = false(nElem,nCases);
            for icase = 1:nCases
                isEdgeCase = edgeCases == obj.code(icase);
                obj.subCellCases(:,icase) = isEdgeCase;
            end
        end     
       
   end
    
   methods (Access = private)
       
       function init(obj,cParams)
           obj.isEdgeCutInElem = cParams.isEdgeCutInElem;
           obj.code   = [5 3 6];  
       end
       
        function d = computeEdgeCases(obj)
            isCut = obj.isEdgeCutInElem;
            nEdges = size(isCut,1);
            edges = (1:nEdges) - 1;
            pow2vector = 2.^(edges);
            d = pow2vector*isCut;
        end
        
   end
        
end