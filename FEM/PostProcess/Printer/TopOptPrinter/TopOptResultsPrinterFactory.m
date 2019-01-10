classdef TopOptResultsPrinterFactory < handle
    
    
    methods (Access = public, Static)
        
        function p = create(printMode,n)
            
            switch printMode
                case 'DesignAndShapes'
                   p = {TopOptShapesPrinter(n),...
                        TopOptDesignVariablePrinter}; 
                case 'DesignVariable'
                   p = {TopOptDesignVariablePrinter};
                case 'DesignAndElementalDensity'
                   p = {TopOptDesignVariablePrinter,...
                       TopOptElementalDensityPrinter}; 
                case 'ElementalDensity'
                   p = {TopOptElementalDensityPrinter};
                case 'DesignElementalDensityAndShape'
                   p = {TopOptDesignVariablePrinter,...
                        TopOptShapesPrinter(n),...
                        TopOptElementalDensityPrinter};
            end                     

        end
        
    end
end