classdef SubCellsCasesComputer < handle
    
   properties (GetAccess = public, SetAccess = private)
       subCellCases
   end
   
   properties (Access = private)
       cutCase
       intergerCodeCases
       integerCases
   end
   
   properties (Access = private)
      connec
      levelSet
   end
   
   methods (Access = public)
       
       function obj = SubCellsCasesComputer(cParams)
          obj.init(cParams);
          obj.computeCutCase();
          obj.computeIntergerCodeCases();
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
           obj.connec  = cParams.connec;
           obj.levelSet = cParams.levelSet;        
       end
       
       function computeIntergerCodeCases(obj)
           switch size(obj.cutCase,2)
               case 3
                  obj.intergerCodeCases = [6 1;5 2;3 4];   
                  %obj.intergerCodeCases = [3 4;5 2;6 1];                        
               case 4
                  obj.intergerCodeCases = [1 14;2 13;4 11;8 7];       
           end           
       end
       
       function computeCutCase(obj)
            nodes = obj.connec;
            ls = zeros(size(nodes));
            for iNode = 1:size(nodes,2)
                ls(:,iNode) = obj.levelSet(nodes(:,iNode));                 
            end            
            obj.cutCase = 1 - heaviside(ls);           
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