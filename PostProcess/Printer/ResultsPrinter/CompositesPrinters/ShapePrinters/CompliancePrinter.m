classdef CompliancePrinter < CompositeResultsPrinter
       
    methods (Access = public)
        
        function obj = CompliancePrinter(d)
            obj.init(d);
        end
        
    end
    
    methods (Access = protected)
        
        function storeFieldsToPrint(obj,d)
            d.fields = d.phyProblems{1}.variables;
            obj.printers{1}.storeFieldsToPrint(d);
            d.fields = d.regDensity;
            obj.printers{2}.storeFieldsToPrint(d);
        end
        
        function createHeadPrinter(obj,d,dh)
            phyPr = d.cost.shapeFunctions{1}.getPhysicalProblems();
            d.quad = phyPr{1}.element.quadrature;
            obj.printers{2}.createHeadPrinter(d,dh);
            h = obj.printers{2}.getHeadPrinter();
            obj.headPrinter = h;            
        end        
        
        function createPrinters(obj,d)
            obj.printers{1} =  ResultsPrinter.create('Elasticity',d);
            obj.printers{2} =  ResultsPrinter.create('DensityGauss',d);
        end       
        
    end
end