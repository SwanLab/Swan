classdef TopOptDesignVariablePrinter < CompositeResultsPrinter
    
    methods (Access = public)
        
        function obj = TopOptDesignVariablePrinter(d)
            obj.init(d);
        end
               
    end
    
    methods (Access = protected)
        
        function createPrinters(obj,d)
            dV = OptimizerToDesignVariableTranslator().translate(d.optimizer);
            obj.printers{1} = ResultsPrinter.create(dV,d);
        end
        
    end
    
end