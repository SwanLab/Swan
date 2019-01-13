classdef TopOptCompliancePrinter < TopOptShapePrinter
    

   methods (Access = public)
              
       function obj = TopOptCompliancePrinter(d,dT)
           obj.simulationStr = d.dStandard.simulationStr;
           obj.createPrinters(d,dT);           
       end
       
   end
   
   methods (Access = private)
       
       function createPrinters(obj,d,dT)
           p =  ResultsPrinter.create('Elasticity',d,dT);
           p.setSimulationStr(obj.simulationStr);
           obj.printers{1} = p;
       end
       
   end
   
end