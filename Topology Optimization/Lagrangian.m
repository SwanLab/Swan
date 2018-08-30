classdef Lagrangian < Objective_Function
    properties
        value_initial 
        lambda
    end
    methods
        function obj=Lagrangian(settings)
            obj.lambda=zeros(1,settings.nconstr);
        end
        function computeFunction(obj,cost,constraint)
%             obj.inactive_constr=-obj.lambda'>constraint.value;
%             constraint.value(obj.inactive_constr)=obj.lambda(obj.inactive_constr);
            obj.value=cost.value + obj.lambda*constraint.value;
        end
        function computeGradient(obj,cost,constraint)
%             constraint.gradient(:,obj.inactive_constr)=0;
            obj.gradient=constraint.gradient*obj.lambda' + cost.gradient;
        end
    end
end
