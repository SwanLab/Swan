classdef ShapePrinter < CompositeResultsPrinter
    
    methods (Access = public)
        
        function obj = ShapePrinter(d)
            obj.init(d);
        end
        
    end
    
    methods (Access = protected)
        
        function storeFieldsToPrint(obj,d)
            for iF = 1:numel(d.fieldToPrint)
                d.fields = d.fieldToPrint{iF}.value;
                obj.printers{iF}.storeFieldsToPrint(d)
            end
        end
        
        function createHeadPrinter(obj,d,dh)
            for iP = 1:numel(obj.printers)
              p = obj.printers{iP};
              hasGauss(iP) = p.getHasGaussData;
            end
            hasGauss = find(hasGauss);
            firstPrinterWithGauss = obj.printers{hasGauss(1)};
            p = firstPrinterWithGauss;
            p.createHeadPrinter(d,dh);
            h = p.getHeadPrinter();
            obj.headPrinter = h;
        end
        
        function createPrinters(obj,d)
            for iF = 1:numel(d.fieldToPrint)
                type = d.fieldToPrint{iF}.type;
                d.name = d.fieldToPrint{iF}.name;
                obj.printers{iF} = ResultsPrinter.create(type,d);
            end
        end
        
    end


end