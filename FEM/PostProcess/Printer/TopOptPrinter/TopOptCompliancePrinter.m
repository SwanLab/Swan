classdef TopOptCompliancePrinter < handle
    
   properties (Access = private)
      simulationStr
      nstre
   end
    
   methods (Access = public)
              
       function obj = TopOptCompliancePrinter(simulationStr)
           obj.simulationStr = simulationStr;
       end
       
       function shPrint = createShapePrinters(obj,d)
           sh =  ElasticityResultsPrinter(d);
           sh.setSimulationStr(obj.simulationStr);
           shPrint{1} = sh;
       end
       
   end
   
end