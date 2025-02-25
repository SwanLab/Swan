classdef Lagrangian < ObjectiveFunction
    
    methods (Access = public)
        
        function obj = Lagrangian(cParams)
            obj.init(cParams);
        end
        
        function computeFunction(obj)
            l  = obj.dualVariable.value;
            c  = obj.constraint.value;
            j  = obj.cost.value;
            obj.value = j + l*c;
        end
        
        function computeGradient(obj)
            l   = obj.dualVariable.value;
            dj  = obj.cost.gradient;
            dc  = obj.constraint.gradient;
            g = dj + dc*l;
            obj.gradient = g;
        end
        
    end
    
    methods (Access = private)
       
    end
    
end
