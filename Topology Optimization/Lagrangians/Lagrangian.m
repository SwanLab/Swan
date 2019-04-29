classdef Lagrangian < ObjectiveFunction
    
    methods (Access = public)
        
        function obj = Lagrangian(cParams)
            obj.init(cParams);
            obj.initializeLambda();
        end
        
        function computeFunction(obj)
            l  = obj.lambda;
            c  = obj.constraint.value;
            j  = obj.cost.value;            
            obj.value = j + l*c;
        end
        
        function computeGradient(obj)
            l   = obj.lambda;
            dj  = obj.cost.gradient;
            dc  = obj.constraint.gradient;
            g = dj + dc*l;
            obj.gradient = g;            
        end
        
    end
    
    methods (Access = private)
       
        function initializeLambda(obj)
            nconstr = numel(obj.constraint);
            obj.lambda = zeros(1,nconstr);            
        end
        
    end
    
end
