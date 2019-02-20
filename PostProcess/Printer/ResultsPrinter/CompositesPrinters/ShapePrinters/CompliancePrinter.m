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
        end
        
        function createPrinters(obj,d)
            p =  ResultsPrinter.create('Elasticity',d);
            obj.printers{1} = p;
        end
        
    end
end