classdef OptimizationMetricsPrinter_AugLag < OptimizationMetricsPrinter
    
    methods (Access = protected)
        
        function printConvergenceVariables(obj,fid)
            fprintf(fid,'Optimality tolerance: %f \n',obj.optimizer.unconstrainedOptimizer.optimalityCond);
            fprintf(fid,'Kappa: %f \n',obj.optimizer.unconstrainedOptimizer.lineSearch.value);
        end
        
    end
    
end

