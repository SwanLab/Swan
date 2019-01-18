classdef TopOptCompliancePrinter < TopOptShapePrinter
    
   methods (Access = public)
              
       function obj = TopOptCompliancePrinter(d)         
           obj.createPrinters(d);           
       end
       
       function itHas = hasGaussData(obj)
           itHas = true;
       end
       
       function storeResultsInfo(obj,compShape)
          phyPr = compShape.getPhysicalProblem();
          d = obj.obtainVariablesAndQuad(phyPr);
          obj.printers{1}.storeResultsInfo(d);
       end
   end
   
   methods (Access = private)
       
       function createPrinters(obj,d)
           p =  ResultsPrinter.create('Elasticity',d);
           obj.printers{1} = p;
       end
       
   end
   
end