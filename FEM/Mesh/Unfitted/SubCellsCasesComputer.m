classdef SubCellsCasesComputer < handle
    
   properties (GetAccess = public, SetAccess = private)
       subCellCases
   end
   
   properties (Access = private)
       cutCase
       intergerCodeCases
       integerCases
   end
   
   methods (Access = public)
       
       function obj = SubCellsCasesComputer(cParams)
          obj.init(cParams);
          obj.computeIntegerCases();                                  
       end
       
       function compute(obj)
            nElem  = size(obj.cutCase,1);
            nCases = size(obj.intergerCodeCases,1);
            obj.subCellCases = false(nElem,nCases);
            for icase = 1:nCases                
                isCaseA = obj.computeIntegerCase(icase,1);
                isCaseB = obj.computeIntegerCase(icase,2);
                isCase = or(isCaseA,isCaseB);
                obj.subCellCases(:,icase) = isCase;
            end
        end     
       
   end
    
   methods (Access = private)
       
       function init(obj,cParams)
           obj.cutCase  = cParams.cutCase;
           switch size(obj.cutCase,2)
               case 3
                  obj.intergerCodeCases = [6 1;5 2;3 4];      
               case 4
                  obj.intergerCodeCases = [1 14;2 13;4 11;8 7];       
           end
       end
       
        function d = computeIntegerCases(obj)
            nodes = obj.cutCase;
            nnode = size(nodes,2);
            nodePos = (1:nnode) - 1;
            pow2vector = 2.^(nodePos);
            d = pow2vector*nodes';
            obj.integerCases = d;
        end
        
        function isCase = computeIntegerCase(obj,icase,subCase)
            integerCase = obj.integerCases;            
            isCase = integerCase == obj.intergerCodeCases(icase,subCase);
        end
        
        
        
   end
        
end