classdef TopOptElementalDensityPrinter < TopOptPrinter
    
    methods (Access = public)
        
        function obj = TopOptElementalDensityPrinter(d)
            p = ResultsPrinter.create('DensityGauss',d);
            obj.printers{1} = p;
        end
        
        function print(obj,istep)
            i = istep;
            obj.printers{1}.printOnlyResults(i);
        end
        
        function storeResultsInfo(obj,d)
            phyPr = d.cost.ShapeFuncs{1}.getPhysicalProblem();
            d.dens = obj.computeElementalDensity(d.x,phyPr);
            d.quad = phyPr.element.quadrature;
            obj.printers{1}.storeResultsInfo(d);
        end
        
        function itHas = hasGaussData(obj)
            itHas = true;
        end                
        
    end
    
    methods (Access = private, Static)
        
        function d = computeElementalDensity(x,phyPr)
            filter = FilterP0(x,phyPr); %Only examples with ls by the moment
            d = filter.getDens0();
        end
        
    end
    
end