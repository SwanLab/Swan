classdef CompositeResultsPrinter < ResultsPrinter
    
    
    properties (Access = protected)
        hasGaussData
        printers
        simulationStr
    end
    
    methods (Access = public)
        
        function  printResults(obj,iter,fileID)
            obj.setSimulationStrToPrinters()
            for iprinter = 1:numel(obj.printers)
                p = obj.printers{iprinter};
                p.printResults(iter,fileID);
            end
        end
        
    end
    
    methods (Access = protected)
        
        function init(obj,d)
            obj.init@ResultsPrinter(d);
            obj.createPrinters(d);
            obj.obtainHasGaussData();
        end
        
        function createHeadPrinter(obj,d,dh)
            for iprinter = 1:numel(obj.printers)
                p = obj.printers{iprinter};
                p.createHeadPrinter(d,dh);
                if p.getHasGaussData()
                    h = p.getHeadPrinter();
                    obj.headPrinter = h;
                    return
                else
                    h = p.getHeadPrinter();
                    obj.headPrinter = h;
                end
            end
        end
        
        function storeFieldsToPrint(obj,d)
            for iprinter = 1:numel(obj.printers)
                p = obj.printers{iprinter};
                p.storeFieldsToPrint(d);
            end
        end
        
    end
    
    methods (Access = private)
        
        function obtainHasGaussData(obj)
            obj.hasGaussData = false;
            for iprinter = 1:numel(obj.printers)
                p = obj.printers{iprinter};
                hG = p.getHasGaussData();
                obj.hasGaussData = hG || obj.hasGaussData;
            end
        end
        
        function setSimulationStrToPrinters(obj)
            s = obj.simulationStr;
            for iprinter = 1:numel(obj.printers)
                p = obj.printers{iprinter};
                p.setSimulationStr(s);
            end
        end
    end
    
    methods (Access = protected, Abstract)
        createPrinters(obj,d)
    end
    
end