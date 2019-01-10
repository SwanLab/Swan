classdef TopOptDesignAndShapesPrinter < TopOptResultsPrinter
       
   properties (Access = protected)
      printers = {TopOptShapesPrinter,TopOptDesignVariablePrinter} 
   end
    
    methods (Access = public)
       
       function obj = TopOptDesignAndShapesPrinter(d)
            obj.compute(d);
       end
   end

end