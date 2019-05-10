classdef OptimizationMetricsPrinter_MMA < OptimizationMetricsPrinter
    
    methods (Access = protected)
        
        function printConvergenceVariables(obj,fid)
            kktnorm = obj.optimizer.historicalVariables.kktnorm;
            fprintf(fid,'Optimality tolerance: %f \n',kktnorm);
        end
        
    end
    
end

