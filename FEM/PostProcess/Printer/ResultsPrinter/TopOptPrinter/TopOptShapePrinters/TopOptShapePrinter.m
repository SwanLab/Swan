classdef TopOptShapePrinter < handle
    
    properties (Access = protected)
        printers
        simulationStr
    end
    
    methods (Access = public)
       
        function setSimulationStr(obj,s)
            for iprinters = 1:numel(obj.printers)
                obj.printers{iprinters}.setSimulationStr(s);
            end
        end
        
        function print(obj,istep)
            for iprinters = 1:numel(obj.printers)
                obj.printers{iprinters}.printOnlyResults(istep);
            end
        end
                
    end
    
    methods (Access = protected, Static)
       
        function d = obtainVariablesAndQuad(phyPr)
            d.quad = phyPr.element.quadrature;
            d.variables = phyPr.variables;                                    
        end
        
    end
    
    methods (Access = public, Abstract)
        hasGaussData(obj)
        storeResultsInfo(obj)
    end
    
end