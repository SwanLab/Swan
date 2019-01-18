classdef TopOptResultsPrinterFactory < handle
    
    properties (Access = private)
        printers
    end
    
    methods (Access = public)
        
        function obj = TopOptResultsPrinterFactory(d)
             obj.buildPrinters(d);
        end
        
        function p = getPrinters(obj)
            p = obj.printers;
        end
        
    end
    
    
     methods (Access = private)
        
        function buildPrinters(obj,d)
            switch d.printMode
                case {'DesignVariable'}                  
                    p = {TopOptDesignVariablePrinter(d)};
                case {'DesignAndElementalDensity'}
                    p = {TopOptDesignVariablePrinter(d),...
                        TopOptElementalDensityPrinter(d)};
                case {'ElementalDensity'}
                    p = {TopOptElementalDensityPrinter(d)};
                case {'DesignAndShapes'}
                    p = {TopOptShapesPrinter(d),...
                        TopOptDesignVariablePrinter(d)};                    
                case {'DesignElementalDensityAndShape'}
                    p = {TopOptDesignVariablePrinter(d),...
                        TopOptShapesPrinter(d),...
                        TopOptElementalDensityPrinter(d)};
            end 
            obj.printers = p;
        end
        
    end
end