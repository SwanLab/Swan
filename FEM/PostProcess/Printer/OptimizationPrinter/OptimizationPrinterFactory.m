classdef OptimizationPrinterFactory < handle
    
    
    methods (Access = public, Static)
        
        function p = create(printMode)
            
            switch printMode
                case {'DesignAndShapes'}
                    p = OptimizerDesignAndShapes();
                case {'DesignVariable'}
                    p = OptimizerPrinterDesignVariable();
                case {'DesignAndElementalDensity'}
                    p = OptimizerDesignAndElementalDensityCase();
                case {'ElementalDensity'}
                    p = OptimizerElemntalDensityCase();
                case {'DesignElementalDensityAndShape'}
                    p = OptimizerDesignElementalDensityAndShape();
            end
            
        end
    end
    
    
end