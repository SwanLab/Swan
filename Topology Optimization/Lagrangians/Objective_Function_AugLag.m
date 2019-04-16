classdef Objective_Function_AugLag < Objective_Function
    
    properties
        value_initial
        lambda
        penalty
    end
    
    properties (Access = private)
        constraintModifier        
    end
    
    methods (Access = public)
        
        function obj=Objective_Function_AugLag(settings)
            obj.lambda=zeros(1,settings.nconstr);
            obj.penalty=ones(1,settings.nconstr);
            obj.createConstraintModifier(settings.constraintCase);
        end
        
        function updateDual(obj,cost,constraint)
            obj.lambda = obj.lambda + obj.penalty.*constraint.value';
            constraint.lambda  = obj.lambda;
            constraint         = obj.modifyInactiveConstraints(constraint);
            obj.computeFunction(cost,constraint);
            obj.computeGradient(cost,constraint);
        end
        
        function updatePrimal(obj,x,cost,constraint)
            cost.computeCostAndGradient(x);
            constraint.computeCostAndGradient(x);
            constraint = obj.modifyInactiveConstraints(constraint);
            obj.computeFunction(cost,constraint)
        end
        
    end
    
    methods (Access = private)
        
        function createConstraintModifier(obj,constraintCase)
            icm = InactiveConstraintsModifier.create(constraintCase);
            obj.constraintModifier = icm;
        end
        
        function computeGradient(obj,cost,constraint)
            l   = obj.lambda;            
            c   = constraint.value;
            dj  = cost.gradient;
            dc  = constraint.gradient;
            rho = obj.penalty;
            g = dj + dc*(l + rho.*c);
            obj.gradient = g;
        end
        
        function computeFunction(obj,cost,constraint)
            l  = obj.lambda;            
            c  = constraint.value;
            j  = cost.value;
            rho = obj.penalty;
            obj.value = j + l*c + 0.5*rho*(c.*c);
        end
        
        function cons = modifyInactiveConstraints(obj,cons)
            obj.constraintModifier.modify(cons,obj.lambda,obj.penalty)
            cons = obj.constraintModifier.constraint;
        end        
        
    end
    
    
end
