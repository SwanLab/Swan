classdef OptimizationMetricsPrinter_Null < OptimizationMetricsPrinter
    
    methods (Access = public)
        
        function print(obj,nIter,iStep)
        end
        
        function printFinal(obj)
        end        
        
    end
    
    methods (Access = protected)
        
        function init(obj,cParams)
        end
        
        function printConvergenceVariables(obj,fid)
        end
        

        
    end
    
end

