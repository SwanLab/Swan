classdef TopOptDesignVariablePrinter < TopOptPrinter
        
    methods (Access = public)         
        
        function create(obj,d)
                  opt = d.optimizer;
            switch opt
                case {'SLERP','PROJECTED SLERP', 'HAMILTON-JACOBI'}
                    p = LevelSetResultsPrinter(d);
                case {'PROJECTED GRADIENT', 'MMA', 'IPOPT'}
                    p = DensityResultsPrinter(d);
            end
            obj.printers = p;
            obj.printers.setSimulationStr(obj.simulationStr);
        end
        
        
        function print(obj,istep,fields)
            i = istep;
            f = fields.designVariable;
            obj.printers.printOnlyResults(i,f);      
        end                               
        
    end
        
    
end