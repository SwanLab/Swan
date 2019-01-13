classdef TopOptShapePrinter < handle
       
   properties (Access = protected)
       printers
       simulationStr      
   end
    
   methods (Access = public, Static)
              
       function p = create(d,dT,shapeName)
            f = TopOptShapePrinterFactory();
            p = f.create(d,dT,shapeName);
       end
       
   end
       
   methods (Access = public)
       
       function p = getPrinters(obj)
            p = obj.printers;
       end      
       
   end
    

end