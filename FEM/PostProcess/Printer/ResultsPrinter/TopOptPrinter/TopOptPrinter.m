classdef TopOptPrinter < handle
    
    properties (Access = protected)
        printers
        simulationStr
    end
    
    methods (Access = public)
                
        function setSimulationStr(obj,s)
            obj.simulationStr = s;
        end
        
    end
    
    methods (Access = public, Abstract)
        print(obj,istep,fields)
    end
end