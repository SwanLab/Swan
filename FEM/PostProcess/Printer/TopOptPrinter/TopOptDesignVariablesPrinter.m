classdef TopOptDesignVariablesPrinter < TopOptResultsPrinter
       
    properties (Access = protected)
        printers = {TopOptDesignVariablePrinter}
    end
    
   methods (Access = public)
       
       function obj = TopOptDesignVariablesPrinter(d)
            obj.compute(d);
       end
   end
   
   methods (Access = protected)
      
%        function createPrinters(obj,d)
%            obj.createDesignVariablePrinter(d);
%        end
%         
%        function  printResults(obj)
%             obj.printDesignVariable();        
%        end
%        
   end
   
end