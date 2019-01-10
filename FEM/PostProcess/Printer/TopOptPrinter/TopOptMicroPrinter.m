classdef TopOptMicroPrinter < handle
    
   properties (Access = private)
      simulationStr
      nstre
   end
    
   methods (Access = public)
              
       function obj = TopOptMicroPrinter(simulationStr,ndim)
           obj.simulationStr = simulationStr;
           if ndim == 2
               obj.nstre = 3;
           elseif ndim == 3
               obj.nstre = 6;
           end
       end
       
       function shPrinter = createShapePrinters(obj,d)
           shPrinter = cell(obj.nstre,1);
           for istre = 1:obj.nstre
               sh = ElasticityMicroResultsPrinter(d);
               sh.setSimulationStr(obj.simulationStr);
               sh.setStrVariablesCase(istre);
               shPrinter{istre} = sh;
           end
                      
       end
       
   end
   
   
    
    
    
    
    
    
    
end