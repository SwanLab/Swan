classdef OptimizationMetricsPrinter_AugLag < OptimizationMetricsPrinter
    
    methods (Access = protected)
        
        function printConvergenceVariables(obj,fid)
            fprintf(fid,'Optimality tolerance: %f \n',obj.optimizer.optimizer_unconstr.opt_cond);
            fprintf(fid,'Kappa: %f \n',obj.optimizer.optimizer_unconstr.line_search.kappa);
        end
        
    end
    
end

