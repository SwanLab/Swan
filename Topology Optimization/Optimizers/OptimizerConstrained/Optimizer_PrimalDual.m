classdef Optimizer_PrimalDual < Optimizer
    
    properties (GetAccess = public, SetAccess = protected)
        unconstrainedOptimizer
    end
    
    properties (Access = protected)
        lagrangian       
        lagrangianSettings
    end
    
    properties (Access = private)

    end
    
    methods (Access = public)
    end
    
    methods (Access = protected)
       
        function updateConvergenceStatus(obj)
            isOptimal   = obj.unconstrainedOptimizer.isOptimal();
            isFeasible  = obj.isFeasible();
            obj.hasConverged = isOptimal && isFeasible;
        end
        
        function itIs = isFeasible(obj)
            active_constr    = true(size(obj.dualVariable.value));
            constraintValues = abs(obj.constraint.value(active_constr));
            constrTol        = obj.targetParameters.constr_tol(active_constr);
            isNotFeasible = any(any(constraintValues > constrTol));
            itIs = ~isNotFeasible;
        end   
        
        function createLagrangian(obj)
            obj.createLagrangianSettings();
            cParams = obj.lagrangianSettings;
            obj.lagrangian = ObjectiveFunction.create(cParams);
        end        
        
       function createOptimizerUnconstrained(obj,cParams)
            cParams.lagrangian      = obj.lagrangian;           
            cParams.convergenceVars = obj.convergenceVars;
            obj.unconstrainedOptimizer = Optimizer_Unconstrained.create(cParams);              
        end             
        
    end
    
    methods (Access = protected, Abstract)
        createLagrangianSettings(obj)        
    end
    
end
