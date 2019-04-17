classdef OptimizationMetricsPrinter_MMA < OptimizationMetricsPrinter
    
    methods (Access = protected)
        
        function printConvergenceVariables(obj,fid)
            fprintf(fid,'Optimality tolerance: %f \n',obj.optimizer.kktnorm);
        end
        
    end
    
end

