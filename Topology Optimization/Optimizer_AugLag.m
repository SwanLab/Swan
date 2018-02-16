classdef Optimizer_AugLag < Optimizer
    properties
        optimizer_unconstr
        penalty
        iter = 0;
    end
    methods
        function obj = Optimizer_AugLag(settings,optimizer_unconstr)
            obj@Optimizer(settings);
            obj.objfunc = Objective_Function_AugLag(settings);
            obj.optimizer_unconstr = optimizer_unconstr;
        end
        function x = updateX(obj,x_ini,cost,constraint,interpolation,filter)
            obj.checkInitial;
            obj.optimizer_unconstr.target_parameters = obj.target_parameters;
            obj.shfunc_volume.target_parameters = obj.target_parameters;
            obj.optimizer_unconstr.shfunc_volume.target_parameters = obj.target_parameters;
            obj.shfunc_volume.computef(x_ini,obj.physicalProblem,interpolation,filter);
            
            obj.objfunc.lambda = obj.objfunc.lambda+obj.objfunc.penalty.*constraint.value';
            obj.objfunc.computeFunction(cost,constraint);
            obj.objfunc.computeGradient(cost,constraint);
            
            obj.optimizer_unconstr.volume_initial = obj.shfunc_volume.value;
            obj.optimizer_unconstr.objfunc = obj.objfunc;
            obj.optimizer_unconstr.objfunc.value_initial = obj.objfunc.value;
            obj.optimizer_unconstr.computeKappa(x_ini,obj.objfunc.gradient);
            obj.optimizer_unconstr.stop_criteria = 1;
            
            obj.optimizer_unconstr.setPhysicalProblem(obj.physicalProblem);
            while (obj.optimizer_unconstr.stop_criteria)
                x = obj.optimizer_unconstr.updateX(x_ini,cost,constraint,interpolation,filter); %x = obj.optimizer_unconstr.updateX(x_ini,cost,constraint,obj.physicalProblem,interpolation,filter);
                obj.stop_vars = obj.optimizer_unconstr.stop_vars;
            end
            
            active_constr = obj.penalty > 0;
            obj.stop_criteria = obj.optimizer_unconstr.opt_cond >=  obj.optimizer_unconstr.optimality_tol || any(any(abs(constraint.value(active_constr)) > obj.optimizer_unconstr.constr_tol(active_constr)));
            
%             % !! For monitoring (NOT DECIDED IF MONITOR AugLag) !!
%             stop_vars{1,1} = obj.optimizer_unconstr.optimality_tol; stop_vars{1,2} = obj.optimizer_unconstr.opt_cond;
%             stop_vars{2,1} = abs(constraint.value(active_constr)); stop_vars{2,2} = obj.optimizer_unconstr.constr_tol(active_constr);
        end
        function checkInitial(obj)
            if isempty(obj.optimizer_unconstr.Ksmooth)
                obj.optimizer_unconstr.Msmooth = obj.Msmooth;
                obj.optimizer_unconstr.Ksmooth = obj.Ksmooth;
                obj.optimizer_unconstr.epsilon_scalar_product_P1 = obj.epsilon_scalar_product_P1;
            end
        end
        
    end
    
end