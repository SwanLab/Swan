classdef SubCellsCasesComputer < handle
    
   properties (GetAccess = public, SetAccess = private)
       subCellCases
       isSubCellsInterior
   end
   
   properties (Access = private)
       isNodeInterior
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
            nElem  = size(obj.isNodeInterior,1);
            nCases = size(obj.intergerCodeCases,1);
            obj.subCellCases = false(nElem,nCases);            
           switch size(obj.isNodeInterior,2)
               case 3
                   nSubCells = 3;
               case 4
                   nSubCells = 4;
           end                                
            obj.isSubCellsInterior = false(nSubCells,nElem);
            for icase = 1:nCases                
                isCaseA = obj.computeIntegerCase(icase,1);
                isCaseB = obj.computeIntegerCase(icase,2);
                isCase = or(isCaseA,isCaseB);
                obj.subCellCases(:,icase) = isCase;
                obj.isSubCellsInterior(2:end,isCaseA) = true;
                obj.isSubCellsInterior(1,isCaseB) = true;
            end
        end     
       
   end
    
   methods (Access = private)
       
       function init(obj,cParams)
           obj.connec  = cParams.connec;
           obj.levelSet = cParams.levelSet;        
       end
       
       function computeIntergerCodeCases(obj)
           switch size(obj.isNodeInterior,2)
               case 3
                  obj.intergerCodeCases = [6 1;5 2;3 4];   
                  %obj.intergerCodeCases = [3 4;5 2;6 1];                        
               case 4
                  obj.intergerCodeCases = [14 1;13 2;11 4;7 8];       
           end           
       end
       
       function computeCutCase(obj)
            nodes = obj.connec;
            ls = zeros(size(nodes));
            for iNode = 1:size(nodes,2)
                ls(:,iNode) = obj.levelSet(nodes(:,iNode));                 
            end            
            obj.isNodeInterior = 1 - heaviside(ls);           
       end
       
        function d = computeIntegerCases(obj)
            nodes = obj.isNodeInterior;
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