classdef L1VectorNormProximal < handle
    
    properties (Access = private)
       lambda 
       imageSize
       designVariable
    end
    
    
    methods (Access = public)
        
        function obj = L1VectorNormProximal(cParams)
            obj.lambda = cParams.lambda;
            obj.designVariable = cParams.designVariable;
            obj.imageSize = cParams.imageSize;
        end
        
        function solve(obj)
            lam = obj.lambda;
            normP = obj.computeNormP();
            no  = max(1,normP/lam);
            p = obj.designVariable.value;
            p = p./[no;no];
            obj.designVariable.value = p;            
        end
        
    end
    
    methods (Access = private)
        
       function normP = computeNormP(obj)
            p  = obj.designVariable.value;
            mn = obj.imageSize.rowsTimesColumns;
            normP = hypot(p(1:mn),p(mn+1:end));
        end 
        
    end
    
    
    
    
    
    
    
    
    
end