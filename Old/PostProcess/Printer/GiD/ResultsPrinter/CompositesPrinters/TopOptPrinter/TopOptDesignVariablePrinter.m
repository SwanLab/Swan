classdef TopOptDesignVariablePrinter < CompositeResultsPrinter
    
    methods (Access = public)
        
        function obj = TopOptDesignVariablePrinter(d)
            obj.init(d);
        end
        
    end
    
    methods (Access = protected)
        
        function createPrinters(obj,d)
            for iV = 1:d.nDesignVariables
                dVi = ['DesignVar',num2str(iV)];
                d.name = dVi;
                type = 'ScalarNodal';
                obj.printers{iV} = ResultsPrinter.create(type,d);
            end
        end
        
        function storeFieldsToPrint(obj,dI)
            for iprinter = 1:numel(obj.printers)
                p = obj.printers{iprinter};
                d.fields = dI.fields{iprinter};
                p.storeFieldsToPrint(d);
            end
        end
        
    end
    
end