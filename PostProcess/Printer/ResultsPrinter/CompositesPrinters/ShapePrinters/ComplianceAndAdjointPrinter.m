classdef ComplianceAndAdjointPrinter < CompositeResultsPrinter
    
    methods (Access = public)
        
        function obj = ComplianceAndAdjointPrinter(d)
            obj.init(d);
        end        
        
    end
    
    methods (Access = protected)
        
        function createPrinters(obj,d)
            obj.printers{1} = obj.createCompliancePrinter(d);
            obj.printers{2} = obj.createAdjointPrinter(d);
        end
        
        function storeFieldsToPrint(obj,d)
            obj.storeFieldsToPrintFromPhyPr(obj.printers{1},d.phyProblems{1});
            obj.storeFieldsToPrintFromPhyPr(obj.printers{2},d.phyProblems{2});
        end
        
    end
    
    methods (Access = private, Static)
        
        function storeFieldsToPrintFromPhyPr(printer,phyPr)
            d.fields = phyPr.variables;            
            printer.storeFieldsToPrint(d);
        end
        
        function p = createCompliancePrinter(d)
            p =  ResultsPrinter.create('Elasticity',d);
            p.setStrVariablesNames('Stress','Strain','Disp');
        end
        
        function p = createAdjointPrinter(d)
            p =  ResultsPrinter.create('Elasticity',d);
            p.setStrVariablesNames('StressAdj','StrainAdj','DispAdj');
        end
    end
end