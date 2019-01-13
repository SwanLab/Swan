classdef TopOptDesignVariablePrinter < TopOptPrinter
            
    methods (Access = public)         
        
        function obj = TopOptDesignVariablePrinter(d,dT,opt)
            dV = designVariable.obtainName(opt);
            obj.printers = ResultsPrinter.create(dV,d,dT);
            obj.printers.setSimulationStr(d.dStandard.simulationStr)
        end
                
        function print(obj,istep,fields)
            i = istep;
            f = fields.designVariable;
            obj.printers.printOnlyResults(i,f);      
        end                               
        
    end
    
end