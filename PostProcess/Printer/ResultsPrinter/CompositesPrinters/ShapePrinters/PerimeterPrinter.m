classdef PerimeterPrinter < CompositeResultsPrinter
       
    methods (Access = public)
        
        function obj = PerimeterPrinter(d)
            obj.simulationStr = 'PerimeterShape';                        
            obj.init(d);
            obj.hasGaussData = false;
        end
        
        
    end
    
    methods (Access = protected)
        
        function storeFieldsToPrint(obj,d)
            d.gradient = d.gradient;
            obj.printers{1}.storeFieldsToPrint(d);
            d.x = d.regDensity;
            obj.printers{2}.storeFieldsToPrint(d);
        end
        
        function createPrinters(obj,d)
            obj.printers{1} =  ResultsPrinter.create('Gradient',d);
            obj.printers{2} =  ResultsPrinter.create('Density',d);
        end       
        
    end
end