classdef TopOptResultsPrinterFactory < handle
    
    properties (Access = private)
        printers
    end
    
    methods (Access = public)
        
        function obj = TopOptResultsPrinterFactory(d,dT,s)
            obj.buildPrinters(d,dT);
            obj.setSimulationStrToPrinters(s);
        end
        
        function p = getPrinters(obj)
            p = obj.printers;
        end
        
    end
    
    
     methods (Access = private)
        
        function buildPrinters(obj,d,dT)
            switch dT.printMode
                case 'DesignAndShapes'
                    s = dT.ShapeNames;
                    opt = dT.optimizer;
                    p = {TopOptShapesPrinter(d,dT,s),...
                        TopOptDesignVariablePrinter(d,dT,opt)};
                case 'DesignVariable'
                    opt = dT.optimizer;                    
                    p = {TopOptDesignVariablePrinter(d,dT,opt)};
                case 'DesignAndElementalDensity'
                    opt = dT.optimizer;                    
                    p = {TopOptDesignVariablePrinter(d,dT,opt),...
                        TopOptElementalDensityPrinter(d,dT)};
                case 'ElementalDensity'
                    p = {TopOptElementalDensityPrinter(d,dT)};
                case 'DesignElementalDensityAndShape'
                    s = dT.ShapeNames;                    
                    opt    = dT.optimizer; 
                    p = {TopOptDesignVariablePrinter(d,dT,opt),...
                        TopOptShapesPrinter(d,dT,s),...
                        TopOptElementalDensityPrinter(d,dT)};
            end 
            obj.printers = p;
        end
        
        function setSimulationStrToPrinters(obj,simulationStr)
            for iprinter = 1:numel(obj.printers)
                p = obj.printers{iprinter};
                p.setSimulationStr(simulationStr);
            end            
        end
    end
end