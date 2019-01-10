classdef TopOptElementalDensityPrinter < TopOptPrinter    
    
    methods (Access = public)

        function create(obj,d)
            p = DensityGaussResultsPrinter(d);
            obj.printers = p;
            obj.printers.setSimulationStr(obj.simulationStr);
        end
        
        function print(obj,istep,fields)
            i = istep;
            dens = fields.dens;
            obj.printers.printOnlyResults(i,dens);     
        end                       
        
    end
        
end