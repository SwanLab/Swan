classdef OptimizationMetricsPrinter_Null < OptimizationMetricsPrinter
    
    methods (Access = public)
        
        function print(obj,nIter,iStep)
        end
        
    end
    
    methods (Access = protected)
        
        function printConvergenceVariables(obj,fid)
        end
        
    end
    
end

