classdef TopOptDesignAndElementalDensityPrinter < TopOptResultsPrinter
    
   properties (Access = protected)
      printers = {TopOptDesignVariablePrinter,TopOptElementalDensityPrinter} 
   end
    
    methods (Access = public)
        
        function obj = TopOptDesignAndElementalDensityPrinter(d)
            obj.compute(d);
        end
    end
    
    methods (Access = protected)
        
       
        
    end
    
end