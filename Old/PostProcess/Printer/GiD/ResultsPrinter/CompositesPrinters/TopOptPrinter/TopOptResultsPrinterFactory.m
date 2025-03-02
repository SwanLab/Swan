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
                    p{1} = TopOptDesignVariablePrinter(d);
                case {'DesignAndElementalDensity'}
                    p{1} = TopOptDesignVariablePrinter(d);
                    p{2} = TopOptElementalDensityPrinter(d);
                case {'ElementalDensity'}
                    p = {TopOptElementalDensityPrinter(d)};
                 case {'DesignAndShapes'}
                     p = {ShapesPrinter(d),...
                         TopOptDesignVariablePrinter(d)};
                case {'DesignElementalDensityAndShape'}
                    p = {TopOptDesignVariablePrinter(d),...
                        ShapesPrinter(d),...
                        TopOptElementalDensityPrinter(d)};
            end 
            obj.printers = p;
        end
        
     end
 
end