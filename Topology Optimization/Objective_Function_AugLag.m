classdef Objective_Function_AugLag < Objective_Function
    properties
        value_initial 
        lambda
        penalty
    end
    methods
        function obj=Objective_Function_AugLag(settings)
            obj.lambda=zeros(1,settings.nconstr);
            obj.penalty=ones(1,settings.nconstr);
        end
        function computeFunction(obj,cost,constraint)
           obj.value=cost.value + obj.lambda*constraint.value + 0.5*obj.penalty*(constraint.value.*constraint.value);
        end
        function computeGradient(obj,cost,constraint)
           obj.gradient=constraint.gradient*obj.lambda' + constraint.gradient*(obj.penalty'.*constraint.value) + cost.gradient;         
        end
    end
end
