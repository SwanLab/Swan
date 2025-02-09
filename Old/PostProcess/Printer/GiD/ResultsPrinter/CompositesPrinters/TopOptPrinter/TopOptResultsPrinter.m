classdef TopOptResultsPrinter < CompositeResultsPrinter
      
    methods (Access = public)
        
        function obj = TopOptResultsPrinter(d)
            obj.simulationStr = 'TopOpt';
            obj.init(d);
        end
        
    end
    
    methods (Access = protected)
        
        function createPrinters(obj,d)
            f = TopOptResultsPrinterFactory(d);
            obj.printers = f.getPrinters();
        end
        
    end
    
end