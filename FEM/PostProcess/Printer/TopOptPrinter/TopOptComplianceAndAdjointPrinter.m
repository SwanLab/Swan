classdef TopOptComplianceAndAdjointPrinter < handle
    
   properties (Access = private)
      simulationStr
      nstre
   end
    
   methods (Access = public)
              
       function obj = TopOptComplianceAndAdjointPrinter(simulationStr)
           obj.simulationStr = simulationStr;
       end
       
       function shPrint = createShapePrinters(obj,d)
           shPrint{1} = obj.createCompliancePrinter(d);
           shPrint{2} = obj.createAdjointPrinter(d);
       end
       
   end
      
   methods (Access = private)
      
       function sh = createCompliancePrinter(obj,d)
           sh =  ElasticityResultsPrinter(d);
           sh.setSimulationStr(obj.simulationStr);
           sh.setStrVariablesNames('Stress','Strain','Disp');
       end
       
       function sh = createAdjointPrinter(obj,d)
           sh = ElasticityResultsPrinter(d);
           sh.setSimulationStr(obj.simulationStr);
           sh.setStrVariablesNames('StressAdj','StrainAdj','DispAdj');
       end
   end
end