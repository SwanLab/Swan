classdef TopOptDesignVariablePrinter < TopOptPrinter
    
    methods (Access = public)
        
        function obj = TopOptDesignVariablePrinter(d)
            dV = designVariable.obtainName(d.optimizer);
            obj.printers{1} = ResultsPrinter.create(dV,d);
        end
        
        function print(obj,istep)
            i = istep;
            obj.printers{1}.printOnlyResults(i);
        end
        
        function itHas =  hasGaussData(obj)
            itHas = false;
        end
        
    end
    
end