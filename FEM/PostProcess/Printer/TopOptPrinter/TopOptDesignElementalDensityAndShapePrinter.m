classdef TopOptDesignElementalDensityAndShapePrinter < TopOptResultsPrinter
    
    properties (Access = protected)
      printers = {TopOptDesignVariablePrinter,...
                  TopOptShapesPrinter,...
                   TopOptElementalDensityPrinter}
    end
    
    methods (Access = public)
        
        function obj = TopOptDesignElementalDensityAndShapePrinter(d)
            obj.compute(d);
        end
    end
    
    
end