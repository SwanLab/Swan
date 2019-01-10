classdef TopOptResultsPrinter < ResultsPrinter
    
    properties (Access = protected)
        simulationStr = 'TopOpt';
    end
    
    properties (Access = private)
        headPrinter
        printers
    end
    
    methods (Access = public)
        
        function obj = TopOptResultsPrinter(d)
            obj.init(d);
            obj.createHeadPrinter();
            obj.printHeader();
            obj.buildPrinters(d.printMode);
            obj.createPrinters(d);
        end
        
    end
    
    methods (Access = protected)
        
        function  printResults(obj)
            for iprinter = 1:numel(obj.printers)
                p = obj.printers{iprinter};
                p.print(obj.istep,obj.fields);
            end
        end
        
        function printHeader(obj)
            
            if obj.hasGaussData
                d.fileID = obj.fileID;
                d.gaussDescriptor = obj.gaussDescriptor;
                d.etype = obj.etype;
                d.ngaus = obj.ngaus;
                d.ndim  = obj.ndim;
                d.posgp = obj.posgp;
                obj.headPrinter.print(d);
            else
                f = obj.fileID;
                obj.headPrinter.print(f);
            end
        end
    end
    
    methods (Access = private)
        
        function createHeadPrinter(obj)
            if obj.hasGaussData
                obj.headPrinter = GaussHeadPrinter;
            else
                obj.headPrinter = NoGaussHeadPrinter;
            end
        end
        
        function buildPrinters(obj,printMode)
            factory = TopOptResultsPrinterFactory();
            p = factory.create(printMode,obj.ndim);
            obj.printers = p;
        end
        
        function createPrinters(obj,d)
            for iprinter = 1:numel(obj.printers)
                p = obj.printers{iprinter};
                p.setSimulationStr(obj.simulationStr);
                p.create(d);
            end
        end
        
    end
    
end