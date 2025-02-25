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
          %  obj.printers{3} = obj.createRegularizedDensityPrinter(d);
        end
        
        function storeFieldsToPrint(obj,d)
            obj.storeComplianceFields(d);
            obj.storeAdjointFields(d);
       %     obj.storeRegularizedDensity(d);
        end

        function createHeadPrinter(obj,d,dh)
            phyPr = d.cost.shapeFunctions{1}.getPhysicalProblems();
            d.quad = phyPr{1}.element.quadrature;
            obj.printers{1}.createHeadPrinter(d,dh);
            h = obj.printers{1}.getHeadPrinter();
            obj.headPrinter = h;
        end
        
    end
    
    methods (Access = private, Static)
        
        function p = createCompliancePrinter(d)
            p =  ResultsPrinter.create('Elasticity',d);
            p.setStrVariablesNames('Stress','Strain','Disp');
        end
        
        function p = createAdjointPrinter(d)
            p =  ResultsPrinter.create('Elasticity',d);
            p.setStrVariablesNames('StressAdj','StrainAdj','DispAdj');
        end
        
     %   function p = createRegularizedDensityPrinter(d)
     %       p = ResultsPrinter.create('DensityGauss',d);
     %   end
        
        function storeFieldsToPrintFromPhyPr(printer,phyPr)
            d.fields = phyPr.variables;
            printer.storeFieldsToPrint(d);
        end
        
    end
    
    methods (Access = private)
        
        function storeComplianceFields(obj,d)
            d.fields = d.phyProblems{1}.variables;
            obj.printers{1}.storeFieldsToPrint(d);
        end
        
        function storeAdjointFields(obj,d)
            d.fields = d.phyProblems{2}.variables;
            obj.printers{2}.storeFieldsToPrint(d);
        end
        
        function storeRegularizedDensity(obj,d)
     %       d.fields = d.regDensity;
     %       obj.printers{3}.storeFieldsToPrint(d);
        end
        
    end
end