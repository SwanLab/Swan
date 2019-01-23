classdef Optimizer_Projected_Slerp < Optimizer_Constrained
    properties
        optimizer_unconstr
        objfunc
        problem
    end
    methods
        function obj = Optimizer_Projected_Slerp(settings,mesh,epsilon)
            obj@Optimizer_Constrained(settings,mesh,settings.monitoring);
            obj.objfunc = Lagrangian(settings);
            obj.optimizer_unconstr = Optimizer_SLERP(settings,epsilon);
        end
        
        function x = updateX(obj,x_ini,cost,constraint)
            x_ini = obj.compute_initial_value(x_ini,cost,constraint);
            obj.updateObjFunc(cost,constraint);
            cost.computeCostAndGradient(x_ini);
            constraint.computeCostAndGradient(x_ini);
            obj.objfunc.computeGradient(cost,constraint);
            obj.objfunc.computeFunction(cost,constraint);
            obj.initUnconstrOpt(x_ini);
            obj.optimizer_unconstr.computeX(x_ini,obj.objfunc.gradient);
            
            obj.has_converged = ~(obj.optimizer_unconstr.opt_cond >=  obj.optimizer_unconstr.optimality_tol);
            if ~obj.has_converged
                x = obj.solveUnconstrainedProblem(x_ini,cost,constraint);
                obj.has_converged = ~(obj.optimizer_unconstr.opt_cond >=  obj.optimizer_unconstr.optimality_tol);
            else
                x = x_ini;
                obj.optimizer_unconstr.stop_vars(1,1) = 0;     obj.optimizer_unconstr.stop_vars(1,2) = obj.optimizer_unconstr.theta;
                obj.optimizer_unconstr.stop_vars(2,1) = 0;     obj.optimizer_unconstr.stop_vars(2,2) = obj.optimizer_unconstr.max_constr_change;
                obj.optimizer_unconstr.stop_vars(3,1) = 1;     obj.optimizer_unconstr.stop_vars(3,2) = obj.optimizer_unconstr.line_search.kappa_min;
                obj.stop_vars = obj.optimizer_unconstr.stop_vars;
            end
        end
        
        function x_ini = compute_initial_value(obj,x_ini,cost,constraint)
            
            obj.problem.solver = 'fzero';
            obj.problem.options = optimset(@fzero);
            obj.problem.objective = @(lambda) obj.compute_feasible_design_variable(lambda,x_ini,cost,constraint,0.1);
            obj.problem.x0 = [0 50];
            obj.problem.options = optimset(obj.problem.options,'TolX',1e-3);
            
            obj.optimizer_unconstr.theta = 0.1;
            lambda = fzero(obj.problem);
            obj.objfunc.lambda = lambda;
            constraint.lambda = obj.objfunc.lambda;
            obj.objfunc.computeGradient(cost,constraint);
            obj.optimizer_unconstr.line_search.initKappa;
            x_ini = obj.optimizer_unconstr.computeX(x_ini,obj.objfunc.gradient);
            
            obj.fhtri= [];
            %obj.plotX(x_ini);
        end
    end
    
    methods (Access = private)
        function x = solveUnconstrainedProblem(obj,x_ini,cost,constraint)
            cost_copy_value = cost.value;
            constraint_copy_value = constraint.value;

            cost_copy_gradient = cost.gradient;
            constraint_copy_gradient = constraint.gradient;
            
            lambda_copy = constraint.lambda;
            
            obj.objfunc.computeGradient(cost,constraint);
            
            obj.optimizer_unconstr.line_search.kfrac = 1.1;
            
            while ~obj.optimizer_unconstr.has_converged
                
                cost.value = cost_copy_value;
                constraint.value = constraint_copy_value;
                
                
                cost.gradient = cost_copy_gradient;
                constraint.gradient = constraint_copy_gradient;
                
                constraint.lambda = lambda_copy;
                obj.objfunc.lambda = lambda_copy;
                
                
                obj.objfunc.computeGradient(cost,constraint);
                obj.objfunc.computeFunction(cost,constraint);
                
                
                obj.problem.objective = @(lambda) obj.compute_feasible_design_variable(lambda,x_ini,cost,constraint,obj.optimizer_unconstr.theta);
                obj.problem.x0 = [0 100];
                lambda = fzero(obj.problem);
                
                obj.objfunc.lambda = lambda;
                constraint.lambda = obj.objfunc.lambda;
                obj.objfunc.computeGradient(cost,constraint);
                x = obj.optimizer_unconstr.computeX(x_ini,obj.objfunc.gradient);
                
                cost.computeCostAndGradient(x);
                obj.objfunc.computeFunction(cost,constraint);
                
                incr_norm_L2  = obj.optimizer_unconstr.norm_L2(x,x_ini);
                incr_cost = (obj.objfunc.value - obj.objfunc.value_initial)/abs(obj.objfunc.value_initial);
                
                obj.optimizer_unconstr.has_converged = ((incr_cost < 0 && incr_norm_L2 < obj.optimizer_unconstr.max_constr_change) || obj.optimizer_unconstr.line_search.kappa <= obj.optimizer_unconstr.line_search.kappa_min);
                
                obj.optimizer_unconstr.stop_vars(1,1) = incr_cost;                                      obj.optimizer_unconstr.stop_vars(1,2) = obj.optimizer_unconstr.theta;
                obj.optimizer_unconstr.stop_vars(2,1) = incr_norm_L2;                                   obj.optimizer_unconstr.stop_vars(2,2) = obj.optimizer_unconstr.max_constr_change;
                obj.optimizer_unconstr.stop_vars(3,1) = obj.optimizer_unconstr.line_search.kappa;       obj.optimizer_unconstr.stop_vars(3,2) = obj.optimizer_unconstr.line_search.kappa_min;
                
                if ~obj.has_converged
                    obj.optimizer_unconstr.line_search.computeKappa;
                end
                
                obj.stop_vars = obj.optimizer_unconstr.stop_vars;
            end
            
            obj.optimizer_unconstr.computeX(x_ini,obj.objfunc.gradient);
        end
        
        function fval = compute_feasible_design_variable(obj,lambda,x_ini,cost,constraint,theta)
            obj.objfunc.lambda = lambda;
            constraint.lambda = obj.objfunc.lambda;
            cost.computeCostAndGradient(x_ini)
            constraint.computeCostAndGradient(x_ini)
            obj.objfunc.computeGradient(cost,constraint);
            x = obj.optimizer_unconstr.computeX(x_ini,obj.objfunc.gradient);
            constraint.computeCostAndGradient(x);
            constraint = obj.setConstraint_case(constraint);
            fval = constraint.value;
        end
        
        function updateObjFunc(obj,cost,constraint)
            obj.optimizer_unconstr.target_parameters = obj.target_parameters;
            obj.objfunc.lambda = obj.objfunc.lambda;
            constraint.lambda = obj.objfunc.lambda;
            constraint =obj.setConstraint_case(constraint);
            obj.objfunc.computeFunction(cost,constraint);
            obj.objfunc.computeGradient(cost,constraint);
        end
        
        function initUnconstrOpt(obj,x_ini)
            obj.optimizer_unconstr.objfunc = obj.objfunc;
            obj.optimizer_unconstr.objfunc.value_initial = obj.objfunc.value;
            obj.optimizer_unconstr.line_search.initKappa(x_ini,obj.objfunc.gradient);
            obj.optimizer_unconstr.has_converged = false;
        end
    end
end