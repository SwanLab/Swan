classdef OptimizationPrinterFactory < handle
    
    
    methods (Access = public, Static)
        
        function p = create(printMode)
                
           switch printMode
                case {'DesignAndShapes','DesignElementalDensityAndShape'}
                     p = OptimizerPrinterWithShapes();
                case {'ElementalDensity','DesignAndElementalDensity'}
                     p = OptimizerPrinterWithGaussData();
                case {'DesignVariable'}
                     p = OptimizerPrinterWithNoGaussData();
           end
           
        end
    end
    
    
end