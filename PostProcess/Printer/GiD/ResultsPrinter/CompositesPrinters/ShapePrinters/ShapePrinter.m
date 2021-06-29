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
            obj.printers{1}.createHeadPrinter(d,dh);
            h = obj.printers{1}.getHeadPrinter();
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