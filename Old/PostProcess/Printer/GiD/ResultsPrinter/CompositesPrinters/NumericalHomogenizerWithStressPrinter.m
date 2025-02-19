classdef NumericalHomogenizerWithStressPrinter < CompositeResultsPrinter
      
    properties (Access = private)
        printerNames
    end
    
    methods (Access = public)
        
        function obj = NumericalHomogenizerWithStressPrinter(d)
            obj.simulationStr = 'NumericalHomogenizerWithStress';
            obj.printerNames = d.printers;
            obj.init(d);
        end
        
    end
    
    methods (Access = protected)
        
        function createPrinters(obj,d)
            for iprint = 1:numel(obj.printerNames)
                printer = obj.printerNames{iprint};
                obj.printers{iprint} = ResultsPrinter.create(printer,d);
            end
        end
        
        function storeFieldsToPrint(obj,d)
            for iprint = 1:numel(obj.printerNames)
                printerName = obj.printerNames{iprint};
                printer  = obj.printers{iprint};
                d.fields = d.var2print{iprint};
                obj.storePrinterFields(printer,printerName,d)
            end
        end
        
    end
    
    methods (Access = private)
        
        function storePrinterFields(obj,printer,printerName,d)
            switch printerName
                case 'NumericalHomogenizer'
                    obj.storeNumericalHomogenizerFields(printer,d);
                case 'HomogenizedTensor'
                    obj.storeMicroFields(printer,d)
            end
        end
    end
    
    methods (Access = private, Static)
        
        function storeNumericalHomogenizerFields(printer,d)
            di.fields = d.fields;
            printer.storeFieldsToPrint(di);
        end
        
        function storeMicroFields(printer,d)
            d.phyProblems{1} = d.fields;
            printer.storeFieldsToPrint(d);
        end
        
    end
    
end