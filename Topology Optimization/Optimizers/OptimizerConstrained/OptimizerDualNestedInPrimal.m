classdef OptimizerDualNestedInPrimal < Optimizer_PrimalDual
    
    properties (GetAccess = public, SetAccess = protected)
        name = 'OptimizerDualNestedInPrimal'
    end
    
    properties (Access = public)
        unconstrainedOptimizer
    end
    
    properties (Access = private)
        desVarChangedValue
        costIncrease
        lagrangian
        constraintProjector
    end
    
    methods (Access = public)
        
        function obj = OptimizerDualNestedInPrimal(cParams)
            obj.init(cParams);
            obj.createLagrangian();
            obj.createOptimizerUnconstrained(cParams.uncOptimizerSettings)
            obj.createConstraintProjector();
        end
        
        
        
        function update(obj)
            obj.updateObjFunc();
            
            obj.unconstrainedOptimizer.init();
            
            obj.cost.updateOld();
            obj.constraint.updateOld();
            obj.designVariable.updateOld();
            obj.dualVariable.updateOld();
            
            
            obj.updateObjFunc();
            obj.lagrangian.updateOld();
            
            
            obj.computeValue();
            
            
            obj.lagrangian.updateOld();
            
            obj.designVariable.updateOld();
            obj.dualVariable.updateOld();
            obj.cost.updateOld()
            obj.constraint.updateOld();
            
            
            while ~obj.hasUnconstraintedOptimizerConverged()
                obj.computeValue();
                obj.unconstrainedOptimizer.line_search.computeKappa();
            end
            obj.updateConvergenceStatus();
        end
        
    end
    
    methods (Access = private)
        
        
        function createLagrangian(obj)
            cParams.cost         = obj.cost;
            cParams.constraint   = obj.constraint;
            cParams.dualVariable = obj.dualVariable;
            obj.lagrangian = Lagrangian(cParams);
        end
        
        function createOptimizerUnconstrained(obj,cParams)
            cParams.lagrangian      = obj.lagrangian;
            cParams.convergenceVars = obj.convergenceVars;
            obj.unconstrainedOptimizer = Optimizer_Unconstrained.create(cParams);
        end
        
        function createConstraintProjector(obj)
            cParams.cost        = obj.cost;
            cParams.constraint  = obj.constraint;
            cParams.designVariable = obj.designVariable;
            cParams.dualVariable = obj.dualVariable;
            cParams.lagrangian  = obj.lagrangian;
            cParams.targetParameters = obj.targetParameters;
            cParams.unconstrainedOptimizer = obj.unconstrainedOptimizer;
            obj.constraintProjector = ConstraintProjector(cParams);
        end
        
        function computeValue(obj)
            obj.constraintProjector.project();
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
        
       
        function updateObjFunc(obj)
            obj.lagrangian.computeFunction();
            obj.lagrangian.computeGradient();
        end
        
        
    end
end