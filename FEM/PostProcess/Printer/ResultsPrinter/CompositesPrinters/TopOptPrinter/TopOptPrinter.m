classdef TopOptPrinter < handle
    
    properties (Access = protected)
        printers
        simulationStr
        fields
    end
    
    methods (Access = public)
        
        function setSimulationStr(obj,s)
            for iprinter = 1:numel(obj.printers)
                p = obj.printers{iprinter};
                p.setSimulationStr(s);
            end
        end
        
        function storeResultsInfo(obj,d)
            for iprinter = 1:numel(obj.printers)
                p = obj.printers{iprinter};
                p.storeResultsInfo(d);
            end
        end
        
    end
    
    methods (Access = public, Abstract)
        print(obj,istep,fields)
    end
    
    methods (Access = public, Abstract)
        hasGaussData(obj)
    end
    
    
end