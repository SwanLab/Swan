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
            obj.unconstrainedOptimizer.init(obj.designVariable.value,obj.objfunc)
        end
        
        function x = update(obj)
            x_ini = obj.designVariable.value;
            x_ini = obj.compute_initial_value(x_ini);
            obj.designVariable.value = x_ini;
            obj.updateObjFunc();
            obj.cost.computeCostAndGradient();
            obj.constraint.computeCostAndGradient();
            obj.objfunc.computeGradient();
            obj.objfunc.computeFunction();
            obj.initUnconstrOpt(x_ini);
            obj.designVariable.value = x_ini;
            obj.unconstrainedOptimizer.compute();
            
            obj.hasConverged = ~(obj.unconstrainedOptimizer.opt_cond >=  obj.unconstrainedOptimizer.optimality_tol);
            if ~obj.hasConverged
                x = obj.solveUnconstrainedProblem(x_ini);
                obj.hasConverged = ~(obj.unconstrainedOptimizer.opt_cond >=  obj.unconstrainedOptimizer.optimality_tol);
            else
                x = x_ini;
                obj.storeConvergedInfo(x_ini,obj.unconstrainedOptimizer)
            end
        end
        
    end
    
    methods (Access = private)
        
        function createLagrangian(obj)
            cParamsL.cost       = obj.cost;
            cParamsL.constraint = obj.constraint;
            obj.objfunc = Lagrangian(cParamsL);           
        end        
        
        function x0 = compute_initial_value(obj,x0)
            
            
            obj.problem.solver = 'fzero';
            obj.problem.options = optimset(@fzero);
            obj.problem.objective = @(lambda) obj.compute_feasible_design_variable(lambda,x0);
            obj.problem.x0 = [0 100];
            obj.problem.options = optimset(obj.problem.options,'TolX',1e-2);
            
            obj.unconstrainedOptimizer.theta = 0.1;
            lambda = fzero(obj.problem);
            obj.objfunc.lambda = lambda;
            obj.constraint.lambda = obj.objfunc.lambda;
            obj.objfunc.computeGradient();
            obj.unconstrainedOptimizer.line_search.initKappa;
            obj.designVariable.value = x0;
            obj.unconstrainedOptimizer.compute();
            
            obj.fhtri = [];
        end
        
        function x = solveUnconstrainedProblem(obj,x0)
            
            obj.costCopy       = obj.cost.clone();
            obj.constraintCopy = obj.constraint.clone();
            
            
            
            obj.unconstrainedOptimizer.line_search.kfrac = 1.1;
            
            while ~obj.unconstrainedOptimizer.hasConverged
                obj.restartCost();
                obj.restartConstraint();
                obj.restartObjFunc();
                obj.designVariable.value = x0;
                
                obj.problem.objective = @(lambda) obj.compute_feasible_design_variable(lambda,x0);
                obj.problem.x0 = [0 1000];
                lambda = fzero(obj.problem);
                
                obj.objfunc.lambda = lambda;
                obj.objfunc.computeGradient();               
                x = obj.unconstrainedOptimizer.compute();
                
                obj.designVariable.value = x;
                obj.cost.computeCostAndGradient();
                obj.objfunc.computeFunction();
                
                incr_norm_L2  = obj.unconstrainedOptimizer.norm_L2(x,x0);
                incr_cost     = obj.objfunc.computeIncrement();
                
                
                obj.desVarChangedValue = incr_norm_L2;
                obj.costIncrease = incr_cost;
                
                obj.storeUnconstrainOptimizerInfo();
                obj.unconstrainedOptimizer.hasConverged = obj.hasUnconstraintedOptimizerConverged();
                
                if ~obj.hasUnconstraintedOptimizerConverged()
                    obj.unconstrainedOptimizer.line_search.computeKappa();
                end
                
                obj.convergenceVars = obj.unconstrainedOptimizer.convergenceVars;
            end
            obj.unconstrainedOptimizer.compute(x0,obj.objfunc.gradient);
        end
        
        function itHas = hasUnconstraintedOptimizerConverged(obj)
            itHas = obj.isStepAcceptable() || obj.isLineSeachTooSmall();
        end
        
        function itIs = isLineSeachTooSmall(obj)
            kappa = obj.unconstrainedOptimizer.line_search.kappa;
            kappa_min = obj.unconstrainedOptimizer.line_search.kappa_min;
            itIs = kappa <= kappa_min;
        end
        
        function itIs = isStepAcceptable(obj)
            costHasDecreased = obj.costIncrease < 0;
            constraintsHavePartiallyChanged = obj.desVarChangedValue < obj.unconstrainedOptimizer.maxIncrNormX;
            itIs = costHasDecreased  && constraintsHavePartiallyChanged;
        end
        
        function storeUnconstrainOptimizerInfo(obj)
            obj.unconstrainedOptimizer.convergenceVars.reset();
            obj.unconstrainedOptimizer.convergenceVars.append(obj.costIncrease);
            obj.unconstrainedOptimizer.convergenceVars.append(obj.desVarChangedValue);
            obj.unconstrainedOptimizer.convergenceVars.append(obj.unconstrainedOptimizer.line_search.kappa);
        end
        
        function fval = compute_feasible_design_variable(obj,lambda,x_ini)
            obj.constraint.lambda = lambda;
            obj.designVariable.value = x_ini;
            obj.cost.computeCostAndGradient();
            obj.constraint.computeCostAndGradient();
            obj.objfunc.lambda    = lambda;            
            obj.objfunc.computeGradient();
            x = obj.unconstrainedOptimizer.compute();
            obj.designVariable.value = x;
            obj.constraint.computeCostAndGradient();
            fval = obj.constraint.value;
        end
        
        function updateObjFunc(obj)
            obj.unconstrainedOptimizer.target_parameters = obj.target_parameters;
            obj.objfunc.computeFunction();
            obj.objfunc.computeGradient();
        end
        
        function initUnconstrOpt(obj,x_ini)
            obj.unconstrainedOptimizer.init(x_ini,obj.objfunc);
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
        end
    end
end