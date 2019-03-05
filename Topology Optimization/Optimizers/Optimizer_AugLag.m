classdef Optimizer_AugLag < Optimizer_Constrained
    
    properties (GetAccess = public, SetAccess = private)
        optimizer_unconstr
        objfunc
        penalty
    end
    
    methods (Access = public)
        
        function obj = Optimizer_AugLag(settings,mesh,optimizer_unconstr)
            obj@Optimizer_Constrained(settings,mesh,settings.monitoring);
            obj.objfunc = Objective_Function_AugLag(settings);
            obj.optimizer_unconstr = optimizer_unconstr;
        end
        
        function x = updateX(obj,x_ini,cost,constraint)
            obj.updateObjFunc(cost,constraint);
            obj.initUnconstrOpt(x_ini);
            
            x = obj.solveUnconstrainedProblem(x_ini,cost,constraint);
            
            obj.updateConvergenceStatus();
        end
        
    end
    
    methods (Access = private)
        
        function x = solveUnconstrainedProblem(obj,x_ini,cost,constraint)
            while ~obj.optimizer_unconstr.has_converged
                x = obj.optimizer_unconstr.updateX(x_ini,cost,constraint);
                obj.stop_vars = obj.optimizer_unconstr.stop_vars;
            end
            
            if ~obj.optimizer_unconstr.designImproved
                x = x_ini;
            end
        end
        
        function updateObjFunc(obj,cost,constraint)
            obj.optimizer_unconstr.target_parameters = obj.target_parameters;
            obj.objfunc.lambda = obj.objfunc.lambda + obj.objfunc.penalty.*constraint.value';
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
        
        function updateConvergenceStatus(obj)
            active_constr = obj.penalty > 0;
            hasNotConverged = obj.optimizer_unconstr.opt_cond >=  obj.optimizer_unconstr.optimality_tol || any(any(abs(constraint.value(active_constr)) > obj.optimizer_unconstr.constr_tol(active_constr)));
            obj.has_converged = ~hasNotConverged;
        end
        
    end
    
end