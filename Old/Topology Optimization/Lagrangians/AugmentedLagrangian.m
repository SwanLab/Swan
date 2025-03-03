classdef AugmentedLagrangian < ObjectiveFunction
    
    properties (Access = public)
        penalty
    end
    
    properties (Access = private)
        constraintModifier
    end
    
    methods (Access = public)
        
        function obj = AugmentedLagrangian(cParams)
            obj.init(cParams);
            obj.createConstraintModifier(cParams);
            obj.createPenalty();
        end
        
        function updateBecauseOfPrimal(obj)
            obj.cost.computeFunctionAndGradient();
            obj.constraint.computeFunctionAndGradient();
            obj.modifyInactiveConstraints();
            obj.computeFunction();
            obj.computeGradient();
        end
        
        function updateBecauseOfDual(obj)
            obj.modifyInactiveConstraints();
            obj.computeFunction();
            obj.computeGradient();
        end

    end
    
    methods (Access = private)
        
        function createPenalty(obj)
            nConstraints = obj.constraint.nSF;
            obj.penalty = ones(1,nConstraints);
        end
        
        function computeFunction(obj)
            l  = obj.dualVariable.value;
            c  = obj.constraint.value;
            j  = obj.cost.value;
            rho = obj.penalty;
            obj.value = j + l*c + 0.5*rho*(c.*c);
        end
        
        function computeGradient(obj)
            l   = obj.dualVariable.value;
            c   = obj.constraint.value;
            dj  = obj.cost.gradient;
            dc  = obj.constraint.gradient;
            rho = obj.penalty;
            g = dj + dc*(l + rho.*c);
            obj.gradient = g;
        end
        
        function modifyInactiveConstraints(obj)
            obj.constraintModifier.modify(obj.penalty)
            obj.constraint = obj.constraintModifier.constraint;
        end
        
        function createConstraintModifier(obj,cParams)
            s.type         = cParams.constraintCase;
            s.dualVariable = cParams.dualVariable;
            s.constraint   = cParams.constraint;
            icm = InactiveConstraintsModifier.create(s);
            obj.constraintModifier = icm;
        end
        
    end
    
end
