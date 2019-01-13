classdef TopOptComplianceAndAdjointPrinter < TopOptShapePrinter
    
   methods (Access = public)
              
       function obj = TopOptComplianceAndAdjointPrinter(d,dT)
           obj.simulationStr = d.dStandard.simulationStr;
           obj.createPrinters(d,dT)
       end
       
   end
      
   methods (Access = private)
       
       function createPrinters(obj,d,dT)
           obj.printers{1} = obj.createCompliancePrinter(d,dT);
           obj.printers{2} = obj.createAdjointPrinter(d,dT);
       end
      
       function p = createCompliancePrinter(obj,d,dT)
           p =  ResultsPrinter.create('Elasticity',d,dT);
           p.setSimulationStr(obj.simulationStr);
           p.setStrVariablesNames('Stress','Strain','Disp');
       end
       
       function p = createAdjointPrinter(obj,d,dT)
           p =  ResultsPrinter.create('Elasticity',d,dT);
           p.setSimulationStr(obj.simulationStr);
           p.setStrVariablesNames('StressAdj','StrainAdj','DispAdj');
       end
   end
end