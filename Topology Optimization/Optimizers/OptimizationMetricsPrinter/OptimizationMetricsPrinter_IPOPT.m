classdef OptimizationMetricsPrinter_IPOPT < OptimizationMetricsPrinter
    
    methods (Access = protected)
        
        function printConvergenceVariables(obj,fid)
            inf_du = obj.optimizer.historicalVariables.inf_du;
            fprintf(fid,'Optimality tolerance: %f \n',inf_du);
        end
        
    end
    
end

