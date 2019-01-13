classdef TopOptResultsPrinter < ResultsPrinter
    
    properties (Access = protected)
        simulationStr = 'TopOpt';
        headPrinter
        headerData  
        printers        
        dResultsPrinter
    end
    
    properties (Access = private)
    end
    
    methods (Access = public, Static)
        
        function p = create(d,dT,hasGaussData)
            f = TopOptResultPrinterFactory();
            p = f.create(d,dT,hasGaussData);            
        end
        
    end
    
    methods (Access = protected)
        
        function init(obj,d,dT)
            obj.init@ResultsPrinter(d.dStandard);
            obj.createHeadPrinter(d);
            obj.createPrinters(d,dT);
        end
        
        function storeDataBase(obj,d)
            fieldsNames = fieldnames(d);
            for ifield = 1:length(fieldsNames)
                fieldName = fieldsNames{ifield};
                fieldValue = d.(fieldName);
                obj.(fieldsNames{ifield}) = fieldValue;
            end
        end
        
        function  printResults(obj)
            for iprinter = 1:numel(obj.printers)
                p = obj.printers{iprinter};
                p.print(obj.iter,obj.fields);
            end
        end
        
        function printHeader(obj)
            obj.headerData.fileID = obj.fileID;
            obj.headPrinter.print(obj.headerData);
        end
        
        function createPrinters(obj,d,dT)
            obj.createDataBaseForResultsPrinter(d)
            dI = obj.dResultsPrinter;
            s = obj.simulationStr;
            factory = TopOptResultsPrinterFactory(dI,dT,s);
            obj.printers = factory.getPrinters();
        end
    end
    
    
    methods (Access = protected, Abstract)
        createDataBaseForResultsPrinter(obj)
        createHeadPrinter(obj)
    end
    
end