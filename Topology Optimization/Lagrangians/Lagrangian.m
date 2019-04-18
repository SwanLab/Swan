classdef Lagrangian < ObjectiveFunction
   
    properties
        lambda
    end
    
    methods (Access = public)
        
        function obj = Lagrangian(settings)
            obj.lambda = zeros(1,settings.nconstr);
        end
        
        function computeFunction(obj,cost,constraint)
            obj.value = cost.value + obj.lambda*constraint.value;
        end
        
        function computeGradient(obj,cost,constraint)
            obj.gradient = constraint.gradient*obj.lambda' + cost.gradient;
        end
        
    end
    
end
