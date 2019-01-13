classdef TopOptElementalDensityPrinter < TopOptPrinter    
    
    methods (Access = public)

        function obj = TopOptElementalDensityPrinter(d,dT)
            p = ResultsPrinter.create('DensityGauss',d,dT);
            p.setSimulationStr(d.dStandard.simulationStr);
            obj.printers = p;
        end
        
        function print(obj,istep,fields)
            i = istep;
            dens = fields.dens;
            obj.printers.printOnlyResults(i,dens);     
        end                       
        
    end
        
end