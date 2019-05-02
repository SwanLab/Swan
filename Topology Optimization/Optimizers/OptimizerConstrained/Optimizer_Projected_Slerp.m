classdef Optimizer_Projected_Slerp < Optimizer_PrimalDual
    
    properties (GetAccess = public, SetAccess = protected)
        name = 'PROJECTED SLERP'
    end
    
    properties (Access = public)
        unconstrainedOptimizer
    end
    
    properties (Access = private)
        desVarChangedValue
        costIncrease
        lagrangian
        problem
        costCopy
        constraintCopy
    end
    
    methods (Access = public)
        
        function obj = Optimizer_Projected_Slerp(cParams)
            obj.init(cParams);
            obj.createLagrangian();
            cParams.uncOptimizerSettings.lagrangian = obj.lagrangian;  
            cParams.uncOptimizerSettings.convergenceVars = obj.convergenceVars;            
            cParams.uncOptimizerSettings.type       = 'SLERP';
            obj.unconstrainedOptimizer = Optimizer_Unconstrained.create(cParams.uncOptimizerSettings);
            obj.unconstrainedOptimizer.line_search.kfrac = 1.05;
            obj.defineProblem();            
        end
        
        function update(obj)
            tolCons = obj.targetParameters.constr_tol;
            obj.problem.options = optimset(obj.problem.options,'TolX',1e-2*tolCons);
           

            obj.unconstrainedOptimizer.init();           
            obj.costCopy                = obj.cost.clone();
            obj.constraintCopy          = obj.constraint.clone();
            obj.designVariable.valueOld = obj.designVariable.value;
            obj.dualVariable.valueOld   = obj.dualVariable.value;
            
            
            obj.updateObjFunc();  
            obj.lagrangian.setInitialValue();                        
            obj.computeValue();
            
            
            obj.lagrangian.setInitialValue();            
            obj.costCopy                = obj.cost.clone();
            obj.constraintCopy          = obj.constraint.clone();
            obj.designVariable.valueOld = obj.designVariable.value;
            obj.dualVariable.valueOld   = obj.dualVariable.value;
     
            
            while ~obj.hasUnconstraintedOptimizerConverged()
                obj.computeValue();
                obj.unconstrainedOptimizer.line_search.computeKappa();
            end
            obj.updateConvergenceStatus();
        end
        
    end
    
    methods (Access = private)
        
        function defineProblem(obj)
            obj.problem.solver = 'fzero';
            obj.problem.options = optimset(@fzero);
            obj.problem.objective = @(lambda) obj.computeFeasibleDesignVariable(lambda);
            obj.problem.x0 = [0 100];
        end
        
        function createLagrangian(obj)
            cParams.cost        = obj.cost;
            cParams.constraint  = obj.constraint;
            cParams.dualVariable = obj.dualVariable;
            obj.lagrangian = Lagrangian(cParams);
        end
        
        function computeValue(obj)
            fzero(obj.problem);
            obj.updateObjFunc();
        end
        
        function updateConvergenceStatus(obj)
            isNotOptimal  = obj.unconstrainedOptimizer.opt_cond >=  obj.unconstrainedOptimizer.optimality_tol;
            isNotFeasible = any(any(abs(obj.constraint.value()) > obj.unconstrainedOptimizer.constr_tol()));
            hasNotConverged = isNotOptimal || isNotFeasible;
            obj.hasConverged = ~hasNotConverged;
        end
        
        function itHas = hasUnconstraintedOptimizerConverged(obj)
            itHas = obj.isStepAcceptable() || obj.isLineSeachTooSmall();
        end
        
        function itIs = isStepAcceptable(obj)
            incr = obj.lagrangian.computeIncrement();
            costHasDecreased = incr < 0;
            itIs = costHasDecreased;
        end
        
        function itIs = isLineSeachTooSmall(obj)
            kappa     = obj.unconstrainedOptimizer.line_search.kappa;
            kappa_min = obj.unconstrainedOptimizer.line_search.kappa_min;
            itIs = kappa <= kappa_min;
        end
               
        function fval = computeFeasibleDesignVariable(obj,lambda)
            obj.restartValues();
            obj.dualVariable.value = lambda;            
            obj.updatePrimalVariableBecauseOfDual();
            obj.constraint.computeCostAndGradient();
            fval = obj.constraint.value;
        end
        
        function updatePrimalVariableBecauseOfDual(obj)
            obj.lagrangian.computeGradient();
            obj.unconstrainedOptimizer.hasConverged = false;            
            obj.unconstrainedOptimizer.compute();
            obj.unconstrainedOptimizer.updateConvergenceParams();
        end
        
        function updateObjFunc(obj)
            obj.lagrangian.computeFunction();
            obj.lagrangian.computeGradient();
        end
        
        function restartValues(obj)
            obj.designVariable.value = obj.designVariable.valueOld;
            obj.dualVariable.value   = obj.dualVariable.valueOld;                        
            obj.restartCost();
            obj.restartConstraint();
            obj.restartObjFunc();            
        end
        
        function restartCost(obj)
            obj.cost.value          = obj.costCopy.value;
            obj.cost.gradient       = obj.costCopy.gradient;
        end
        
        function restartConstraint(obj)
            obj.constraint.value    = obj.constraintCopy.value;
            obj.constraint.gradient = obj.constraintCopy.gradient;
        end
        
        function restartObjFunc(obj)
            obj.lagrangian.computeGradient();
            obj.lagrangian.computeFunction();
        end
    end
end