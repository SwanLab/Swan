classdef TopOptMicroPrinter < TopOptShapePrinter
    
   properties (Access = private)
      nstre
   end
    
   methods (Access = public)
              
       function obj = TopOptMicroPrinter(d,dT)
           obj.simulationStr = d.dStandard.simulationStr;           
           obj.computeNstre(d.dStandard.ndim);
           obj.createPrinters(d,dT);
       end
   
   end
   
   methods (Access = private)
       
       function computeNstre(obj,ndim)
           if ndim == 2
               obj.nstre = 3;
           elseif ndim == 3
               obj.nstre = 6;
           end
       end
       
       function createPrinters(obj,d,dT)
           obj.printers = cell(obj.nstre,1);
           for istre = 1:obj.nstre
               sh = ResultsPrinter.create('ElasticityMicro',d,dT);
               sh.setSimulationStr(obj.simulationStr);
               sh.setStrVariablesCase(istre);
               obj.printers{istre} = sh;
           end                      
       end
       
   end

end