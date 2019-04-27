classdef Optimizer_Projected_Slerp < Optimizer_Constrained
    
    properties (GetAccess = public, SetAccess = protected)
        name = 'PROJECTED SLERP'
    end
    
    properties
        unconstrainedOptimizer
        objfunc
        problem
    end
    
    properties (Access = private)
        desVarChangedValue
        costIncrease
    end
    
    methods (Access = public)
        
        function obj = Optimizer_Projected_Slerp(cParams)
            obj.init(cParams);
            obj.objfunc = Lagrangian(cParams);
            obj.unconstrainedOptimizer = Optimizer_SLERP(cParams.uncOptimizerSettings);
        end
        
        function x = update(obj)
            x_ini = obj.designVar.value;
            cost       = obj.cost;
            constraint = obj.constraint;
            x_ini = obj.compute_initial_value(x_ini);
            obj.updateObjFunc();
            cost.computeCostAndGradient(x_ini);
            constraint.computeCostAndGradient(x_ini);
            obj.objfunc.computeGradient(cost,constraint);
            obj.objfunc.computeFunction(cost,constraint);
            obj.initUnconstrOpt(x_ini);
            obj.unconstrainedOptimizer.compute(x_ini,obj.objfunc.gradient);
            
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
        
        function x0 = compute_initial_value(obj,x0)
            
            cost       = obj.cost;
            constraint = obj.constraint;
            
            obj.problem.solver = 'fzero';
            obj.problem.options = optimset(@fzero);
            obj.problem.objective = @(lambda) obj.compute_feasible_design_variable(lambda,x0,cost,constraint);
            obj.problem.x0 = [0 100];
            obj.problem.options = optimset(obj.problem.options,'TolX',1e-2);
            
            obj.unconstrainedOptimizer.theta = 0.1;
            lambda = fzero(obj.problem);
            obj.objfunc.lambda = lambda;
            constraint.lambda = obj.objfunc.lambda;
            obj.objfunc.computeGradient(obj.cost,obj.constraint);
            obj.unconstrainedOptimizer.line_search.initKappa;
            x0 = obj.unconstrainedOptimizer.compute(x0,obj.objfunc.gradient);
            
            obj.fhtri = [];
        end
        
        
        
        function x = solveUnconstrainedProblem(obj,x0)
            cost       = obj.cost;
            constraint = obj.constraint;
            cost_copy_value       = cost.value;
            constraint_copy_value = constraint.value;
            
            cost_copy_gradient       = cost.gradient;
            constraint_copy_gradient = constraint.gradient;
            
            lambda_copy = constraint.lambda;
            
            obj.objfunc.computeGradient(cost,constraint);
            
            obj.unconstrainedOptimizer.line_search.kfrac = 1.1;
            
            while ~obj.unconstrainedOptimizer.hasConverged
                
                cost.value = cost_copy_value;
                constraint.value = constraint_copy_value;
                
                cost.gradient       = cost_copy_gradient;
                constraint.gradient = constraint_copy_gradient;
                
                constraint.lambda   = lambda_copy;
                obj.objfunc.lambda  = lambda_copy;
                
                
                obj.objfunc.computeGradient(cost,constraint);
                obj.objfunc.computeFunction(cost,constraint);
                
                
                obj.problem.objective = @(lambda) obj.compute_feasible_design_variable(lambda,x0,cost,constraint);
                obj.problem.x0 = [0 1000];
                lambda = fzero(obj.problem);
                
                obj.objfunc.lambda = lambda;
                obj.objfunc.computeGradient(cost,constraint);
                x = obj.unconstrainedOptimizer.compute(x0,obj.objfunc.gradient);
                
                cost.computeCostAndGradient(x);
                obj.objfunc.computeFunction(cost,constraint);
                
                incr_norm_L2  = obj.unconstrainedOptimizer.norm_L2(x,x0);
                incr_cost = obj.objfunc.computeIncrement();
                
                
                obj.desVarChangedValue = incr_norm_L2;
                obj.costIncrease = incr_cost;
                
                obj.storeUnconstrainOptimizerInfo();
                obj.unconstrainedOptimizer.hasConverged = obj.hasUnconstraintedOptimizerConverged();
                
                if ~obj.hasUnconstraintedOptimizerConverged()
                    obj.unconstrainedOptimizer.line_search.computeKappa;
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
        
        function fval = compute_feasible_design_variable(obj,lambda,x_ini,cost,constraint)
            obj.objfunc.lambda = lambda;
            constraint.lambda = obj.objfunc.lambda;
            cost.computeCostAndGradient(x_ini)
            constraint.computeCostAndGradient(x_ini)
            obj.objfunc.computeGradient(cost,constraint);
            x = obj.unconstrainedOptimizer.compute(x_ini,obj.objfunc.gradient);
            constraint.computeCostAndGradient(x);
            fval = constraint.value;
        end
        
        function updateObjFunc(obj)
            obj.unconstrainedOptimizer.target_parameters = obj.target_parameters;
            obj.objfunc.lambda = obj.objfunc.lambda;
            obj.constraint.lambda = obj.objfunc.lambda;
            obj.objfunc.computeFunction(obj.cost,obj.constraint);
            obj.objfunc.computeGradient(obj.cost,obj.constraint);
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
        
    end
end