classdef InactiveConstraintsModifier < handle
    
    properties (Access = public)
        constraint
    end
    
    properties (Access = private)
        dualVariable
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
        
        function obj = InactiveConstraintsModifier(cParams)
           obj.dualVariable = cParams.dualVariable;
           obj.constraint   = cParams.constraint;
        end
        
       function modify(obj,penalty)
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
            t = -obj.dualVariable.value(:)./obj.penalty(:);
            obj.threshold = t';
        end
        
        function modifyValue(obj)
            obj.constraint.value(obj.isInactive) = obj.threshold;
        end
        
        function modifyGradient(obj)
            obj.constraint.gradient(:,obj.isInactive) = 0;
        end
        
    end
    
end