classdef Optimizer_Projected_Slerp < Optimizer_Constrained
    
    properties
        optimizer_unconstr
        objfunc
        problem
    end
    
    properties (Access = private)
        desVarChangedValue
        costIncrease
    end
    
    methods (Access = public)
        
        function obj = Optimizer_Projected_Slerp(settings,mesh)
            
            ocS.settings        = settings;
            ocS.designVariable  = settings.designVar;
            ocS.monitoring      = settings.monitoring;
            
            obj@Optimizer_Constrained(ocS);%settings,mesh,settings.monitoring);
            obj.objfunc = Lagrangian(settings);
            obj.optimizer_unconstr = Optimizer_SLERP(settings.uncOptimizerSettings);
        end
        
        function x = update(obj,x_ini)
            cost       = obj.cost;
            constraint = obj.constraint;
            x_ini = obj.compute_initial_value(x_ini);
            obj.updateObjFunc();
            cost.computeCostAndGradient(x_ini);
            constraint.computeCostAndGradient(x_ini);
            obj.objfunc.computeGradient(cost,constraint);
            obj.objfunc.computeFunction(cost,constraint);
            obj.initUnconstrOpt(x_ini);
            obj.optimizer_unconstr.computeX(x_ini,obj.objfunc.gradient);
            
            obj.hasConverged = ~(obj.optimizer_unconstr.opt_cond >=  obj.optimizer_unconstr.optimality_tol);
            if ~obj.hasConverged
                x = obj.solveUnconstrainedProblem(x_ini);
                obj.hasConverged = ~(obj.optimizer_unconstr.opt_cond >=  obj.optimizer_unconstr.optimality_tol);
            else
                x = x_ini;
                obj.storeConvergedInfo(x_ini,obj.optimizer_unconstr)
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
            
            obj.optimizer_unconstr.theta = 0.1;
            lambda = fzero(obj.problem);
            obj.objfunc.lambda = lambda;
            constraint.lambda = obj.objfunc.lambda;
            obj.objfunc.computeGradient(obj.cost,obj.constraint);
            obj.optimizer_unconstr.line_search.initKappa;
            x0 = obj.optimizer_unconstr.computeX(x0,obj.objfunc.gradient);
            
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
            
            obj.optimizer_unconstr.line_search.kfrac = 1.1;
            
            while ~obj.optimizer_unconstr.hasConverged
                
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
                x = obj.optimizer_unconstr.computeX(x0,obj.objfunc.gradient);
                
                cost.computeCostAndGradient(x);
                obj.objfunc.computeFunction(cost,constraint);
                
                incr_norm_L2  = obj.optimizer_unconstr.norm_L2(x,x0);
                incr_cost = obj.objfunc.computeIncrement();
                
                
                obj.desVarChangedValue = incr_norm_L2;
                obj.costIncrease = incr_cost;
                
                obj.storeUnconstrainOptimizerInfo();
                obj.optimizer_unconstr.hasConverged = obj.hasUnconstraintedOptimizerConverged();
                
                if ~obj.hasUnconstraintedOptimizerConverged()
                    obj.optimizer_unconstr.line_search.computeKappa;
                end
                
                obj.stop_vars = obj.optimizer_unconstr.stop_vars;
            end
            
            obj.optimizer_unconstr.computeX(x0,obj.objfunc.gradient);
        end
        
        function itHas = hasUnconstraintedOptimizerConverged(obj)
            itHas = obj.isStepAcceptable() || obj.isLineSeachTooSmall();
        end
        
        function itIs = isLineSeachTooSmall(obj)
            kappa = obj.optimizer_unconstr.line_search.kappa;
            kappa_min = obj.optimizer_unconstr.line_search.kappa_min;
            itIs = kappa <= kappa_min;
        end
        
        function itIs = isStepAcceptable(obj)
            costHasDecreased = obj.costIncrease < 0;
            constraintsHavePartiallyChanged = obj.desVarChangedValue < obj.optimizer_unconstr.maxIncrNormX;
            itIs = costHasDecreased  && constraintsHavePartiallyChanged;
        end
        
        function storeUnconstrainOptimizerInfo(obj)
            obj.optimizer_unconstr.stop_vars(1,1) = obj.costIncrease;
            obj.optimizer_unconstr.stop_vars(1,2) = obj.optimizer_unconstr.theta;
            obj.optimizer_unconstr.stop_vars(2,1) = obj.desVarChangedValue;
            obj.optimizer_unconstr.stop_vars(2,2) = obj.optimizer_unconstr.maxIncrNormX;
            obj.optimizer_unconstr.stop_vars(3,1) = obj.optimizer_unconstr.line_search.kappa;
            obj.optimizer_unconstr.stop_vars(3,2) = obj.optimizer_unconstr.line_search.kappa_min;
        end
        
        function fval = compute_feasible_design_variable(obj,lambda,x_ini,cost,constraint)
            obj.objfunc.lambda = lambda;
            constraint.lambda = obj.objfunc.lambda;
            cost.computeCostAndGradient(x_ini)
            constraint.computeCostAndGradient(x_ini)
            obj.objfunc.computeGradient(cost,constraint);
            x = obj.optimizer_unconstr.computeX(x_ini,obj.objfunc.gradient);
            constraint.computeCostAndGradient(x);
            % constraint = obj.objfunc.setConstraint_case(constraint,obj.constraint_case);
            fval = constraint.value;
        end
        
        function updateObjFunc(obj)
            obj.optimizer_unconstr.target_parameters = obj.target_parameters;
            obj.objfunc.lambda = obj.objfunc.lambda;
            obj.constraint.lambda = obj.objfunc.lambda;
            % constraint =obj.setConstraint_case(constraint);
            obj.objfunc.computeFunction(obj.cost,obj.constraint);
            obj.objfunc.computeGradient(obj.cost,obj.constraint);
        end
        
        function initUnconstrOpt(obj,x_ini)
            obj.optimizer_unconstr.init(x_ini,obj.objfunc);
        end
        
        %         function initUnconstrOpt(obj,x_ini)
        %             obj.optimizer_unconstr.objfunc = obj.objfunc;
        %             obj.optimizer_unconstr.objfunc.value_initial = obj.objfunc.value;
        %             obj.optimizer_unconstr.line_search.initKappa(x_ini,obj.objfunc.gradient);
        %             obj.optimizer_unconstr.hasConverged = false;
        %         end
        %
        function storeConvergedInfo(obj)
            obj.optimizer_unconstr.stop_vars(1,1) = 0;     obj.optimizer_unconstr.stop_vars(1,2) = obj.optimizer_unconstr.theta;
            obj.optimizer_unconstr.stop_vars(2,1) = 0;     obj.optimizer_unconstr.stop_vars(2,2) = obj.optimizer_unconstr.maxIncrNormX;
            obj.optimizer_unconstr.stop_vars(3,1) = 1;     obj.optimizer_unconstr.stop_vars(3,2) = obj.optimizer_unconstr.line_search.kappa_min;
            obj.stop_vars = obj.optimizer_unconstr.stop_vars;
        end
        
    end
end