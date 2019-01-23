classdef TopOptDesignVariablePrinter < CompositeResultsPrinter
    
    methods (Access = public)
        
        function obj = TopOptDesignVariablePrinter(d)
            obj.init(d);
        end
               
    end
    
    methods (Access = protected)
        
        function createPrinters(obj,d)
            dV = designVariable.obtainName(d.optimizer);
            obj.printers{1} = ResultsPrinter.create(dV,d);
        end
        
%         function storeFieldsToPrint(obj,d)
%             obj.printers{1}.storeFieldsToPrint(d)            
%         end        
        
    end
    
end