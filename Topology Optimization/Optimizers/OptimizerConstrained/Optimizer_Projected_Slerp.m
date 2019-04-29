classdef Optimizer_Projected_Slerp < Optimizer_Constrained
    
    properties (GetAccess = public, SetAccess = protected)
        name = 'PROJECTED SLERP'
    end
    
    properties (Access = public)
        unconstrainedOptimizer
    end
    
    properties (Access = private)
        desVarChangedValue
        costIncrease
        objfunc
        problem
        costCopy
        constraintCopy
    end
    
    methods (Access = public)
        
        function obj = Optimizer_Projected_Slerp(cParams)
            obj.init(cParams);
            obj.createLagrangian();
            obj.unconstrainedOptimizer = Optimizer_SLERP(cParams.uncOptimizerSettings);
            obj.unconstrainedOptimizer.init(obj.objfunc)
            obj.defineProblem();
        end
        
        function update(obj)
            obj.unconstrainedOptimizer.line_search.initKappa;
            obj.unconstrainedOptimizer.line_search.kfrac = 1.1;
            
            nConstraints = obj.constraint.nSF;
            obj.constraint.lambda  = zeros(1,nConstraints);
            obj.costCopy           = obj.cost.clone();
            obj.constraintCopy     = obj.constraint.clone();
            obj.designVariable.valueOld = obj.designVariable.value;

            obj.computeValue();
            
            
            obj.objfunc.setInitialValue();
            
            while ~obj.hasUnconstraintedOptimizerConverged()
                obj.unconstrainedOptimizer.line_search.computeKappa();
                obj.computeValue();
            end
            obj.updateConvergenceStatus()
        end
        
    end
    
    methods (Access = private)
        
        function defineProblem(obj)
            obj.problem.solver = 'fzero';
            obj.problem.options = optimset(@fzero);
            obj.problem.objective = @(lambda) obj.computeFeasibleDesignVariable(lambda);
            obj.problem.x0 = [0 100];
            obj.problem.options = optimset(obj.problem.options,'TolX',1e-4);
        end
        
        function createLagrangian(obj)
            cParamsL.cost       = obj.cost;
            cParamsL.constraint = obj.constraint;
            obj.objfunc = Lagrangian(cParamsL);
        end
        
        function computeValue(obj)
            obj.unconstrainedOptimizer.theta = 0.1;
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
            incr = obj.objfunc.computeIncrement();
            costHasDecreased = incr < 0;
            constraintsHavePartiallyChanged = true;
            %constraintsHavePartiallyChanged = obj.desVarChangedValue < obj.unconstrainedOptimizer.maxIncrNormX;
            itIs = costHasDecreased  && constraintsHavePartiallyChanged;
        end
        
        function itIs = isLineSeachTooSmall(obj)
            kappa     = obj.unconstrainedOptimizer.line_search.kappa;
            kappa_min = obj.unconstrainedOptimizer.line_search.kappa_min;
            itIs = kappa <= kappa_min;
        end
        
        function storeUnconstrainOptimizerInfo(obj)
            obj.unconstrainedOptimizer.convergenceVars.reset();
            obj.unconstrainedOptimizer.convergenceVars.append(obj.costIncrease);
            obj.unconstrainedOptimizer.convergenceVars.append(obj.desVarChangedValue);
            obj.unconstrainedOptimizer.convergenceVars.append(obj.unconstrainedOptimizer.line_search.kappa);
        end
        
        function fval = computeFeasibleDesignVariable(obj,lambda)
            obj.designVariable.value = obj.designVariable.valueOld;
            obj.restartCost();
            obj.restartConstraint();
            obj.restartObjFunc();
            obj.updatePrimalVariableBecauseOfDual(lambda);
            obj.constraint.computeCostAndGradient();
            fval = obj.constraint.value;
        end
        
        function updatePrimalVariableBecauseOfDual(obj,lambda)
            obj.constraint.lambda = lambda;
            obj.unconstrainedOptimizer.init3();
            obj.objfunc.lambda    = lambda;
            obj.objfunc.computeGradient();
            obj.unconstrainedOptimizer.compute();
        end
        
        function updateObjFunc(obj)
            obj.unconstrainedOptimizer.target_parameters = obj.target_parameters;
            obj.objfunc.computeFunction();
            obj.objfunc.computeGradient();
        end
        
        
        function storeConvergedInfo(obj)
            obj.unconstrainedOptimizer.convergenceVars.reset();
            obj.unconstrainedOptimizer.convergenceVars.append(obj.costIncrease);
            obj.unconstrainedOptimizer.convergenceVars.append(obj.desVarChangedValue);
            obj.unconstrainedOptimizer.convergenceVars.append(obj.unconstrainedOptimizer.line_search.kappa);
            obj.convergenceVars = obj.unconstrainedOptimizer.convergenceVars;
        end
        
        function restartCost(obj)
            obj.cost.value          = obj.costCopy.value;
            obj.cost.gradient       = obj.costCopy.gradient;
        end
        
        function restartConstraint(obj)
            obj.constraint.value    = obj.constraintCopy.value;
            obj.constraint.gradient = obj.constraintCopy.gradient;
            obj.constraint.lambda   = obj.constraintCopy.lambda;
        end
        
        function restartObjFunc(obj)
            obj.objfunc.lambda      = obj.constraintCopy.lambda;
            obj.objfunc.computeGradient();
            obj.objfunc.computeFunction();
            obj.objfunc.setInitialValue();            
        end
    end
end