classdef TopOptElementalDensitysPrinter < TopOptResultsPrinter
    
    properties (Access = protected)
        printers = {TopOptElementalDensityPrinter}
    end
    
    
    methods (Access = public)
        
        function obj = TopOptElementalDensitysPrinter(d)
            obj.compute(d);
        end
    end
    
    
end