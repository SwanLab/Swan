classdef InactiveConstraintsModifier < handle
    
    properties (Access = public)
        constraint        
    end
    
    properties (Access = private)
        lambda
        penalty
        isInactive
        threshold        
    end    
    
    methods (Access = public, Static)
        
        function obj = create(constraintCase)
            factory = InactiveConstraintsModifierFactory;
            obj = factory.create(constraintCase);            
        end
        
    end
    
    methods (Access = public)
        
       function modify(obj,cons,lambda,penalty)
            obj.constraint = cons;
            obj.lambda = lambda;
            obj.penalty = penalty;            
            obj.computeThreshold();
            obj.obtainInactiveConstraints();
            obj.modifyValue();
            obj.modifyGradient();            
        end
        
    end
    
    methods (Access = private)
                               
        function obtainInactiveConstraints(obj)
            obj.isInactive = obj.threshold > obj.constraint.value;
        end
        
        function computeThreshold(obj)
            t = -obj.lambda(:)./obj.penalty(:);
            obj.threshold = t';
        end
        
        function modifyValue(obj)
            obj.constraint.value(obj.isInactive);
        end
        
        function modifyGradient(obj)
            obj.constraint.gradient(:,obj.isInactive) = 0;
        end
        
    end
    
end