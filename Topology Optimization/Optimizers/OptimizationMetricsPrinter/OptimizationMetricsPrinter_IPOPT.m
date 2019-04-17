classdef OptimizationMetricsPrinter_IPOPT < OptimizationMetricsPrinter
    
    methods (Access = protected)
        
        function printConvergenceVariables(obj,fid)
            fprintf(fid,'Optimality tolerance: %f \n',obj.optimizer.data.inf_du);
        end
        
    end
    
end

