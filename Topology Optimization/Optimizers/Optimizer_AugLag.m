classdef Optimizer_AugLag < Optimizer_Constrained
    properties
        optimizer_unconstr
        objfunc
        penalty
    end
    methods
        function obj = Optimizer_AugLag(settings,mesh,optimizer_unconstr)
            obj@Optimizer_Constrained(settings,mesh,settings.monitoring);
            obj.objfunc = Objective_Function_AugLag(settings);
            obj.optimizer_unconstr = optimizer_unconstr;
        end
        
        function x = updateX(obj,x_ini,cost,constraint)
            obj.updateObjFunc(cost,constraint);
            obj.initUnconstrOpt(x_ini);
            
            x = obj.solveUnconstrainedProblem(x_ini,cost,constraint);
            
            active_constr = obj.penalty > 0;
            obj.stop_criteria = obj.optimizer_unconstr.opt_cond >=  obj.optimizer_unconstr.optimality_tol || any(any(abs(constraint.value(active_constr)) > obj.optimizer_unconstr.constr_tol(active_constr)));
        end
    end
    
    methods (Access = private)
        function x = solveUnconstrainedProblem(obj,x_ini,cost,constraint)
            while obj.optimizer_unconstr.stop_criteria
                x = obj.optimizer_unconstr.updateX(x_ini,cost,constraint); %x = obj.optimizer_unconstr.updateX(x_ini,cost,constraint,obj.physicalProblem,interpolation,filter);
                obj.stop_vars = obj.optimizer_unconstr.stop_vars;
            end
        end
        
        function updateObjFunc(obj,cost,constraint)
            obj.optimizer_unconstr.target_parameters = obj.target_parameters;
            obj.objfunc.lambda = obj.objfunc.lambda + obj.objfunc.penalty.*constraint.value';
            constraint.lambda = obj.objfunc.lambda;
            obj.objfunc.computeFunction(cost,constraint);
            obj.objfunc.computeGradient(cost,constraint);
        end
        
        function initUnconstrOpt(obj,x_ini)
            obj.optimizer_unconstr.objfunc = obj.objfunc;
            obj.optimizer_unconstr.objfunc.value_initial = obj.objfunc.value;
            obj.optimizer_unconstr.computeKappa(x_ini,obj.objfunc.gradient);
            obj.optimizer_unconstr.stop_criteria = 1;
        end
    end
end