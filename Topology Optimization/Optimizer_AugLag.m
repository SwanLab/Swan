classdef Optimizer_AugLag < Optimizer
    properties
        lambda
        max_constr_change
        stop_Criteria_x
        kappa_min
        shfunc_volume
        optimizer_unconstr
        penalty  
        iter=0;
    end 
    methods
        function obj=Optimizer_AugLag(settings,optimizer_unconstr)
            obj@Optimizer(settings);
            obj.lambda=0;
            obj.optimizer_unconstr=optimizer_unconstr;
            obj.penalty=ones(settings.nconstr,1);
            obj.max_constr_change=+Inf;
            obj.shfunc_volume=ShFunc_Volume(settings.volume);
            obj.kappa_min=1.0000e-15;
        end
        function x_ls=updateX(obj,x_ini,cost,constraint, physProblem, interpolation,filter)  
            obj.checkInitial;
            obj.iter=obj.iter+1;
            iter=obj.iter
            obj.shfunc_volume.computef(x_ini,physProblem,interpolation,filter);
            volume = obj.shfunc_volume.value;
            obj.lambda = obj.lambda+obj.penalty*constraint.value;
            cost_ini = cost.value + obj.lambda*constraint.value + 0.5*obj.penalty*(constraint.value.*constraint.value);
            gradient_ini = constraint.gradient*obj.lambda' + constraint.gradient*(obj.penalty'.*constraint.value) + cost.gradient;
           
            obj.optimizer_unconstr.computeKappa(x_ini,gradient_ini);
            obj.stop_Criteria_x=1;          
            while(obj.stop_Criteria_x)
                x_ls=obj.optimizer_unconstr.updateX(x_ini,gradient_ini);
                physProblem=obj.updateEquilibrium(x_ls,physProblem,interpolation,filter);
                cost.computef(x_ls,physProblem,interpolation,filter);
                constraint.computef(x_ls,physProblem,interpolation,filter);
                obj.shfunc_volume.computef(x_ls,physProblem,interpolation,filter);
                
                cost_ls = cost.value + obj.lambda*constraint.value + 0.5*obj.penalty*(constraint.value.*constraint.value);
                volume_ls =obj.shfunc_volume.value;
                
                incr_vol_ls = abs(volume_ls - volume);
                incr_cost = (cost_ls - cost_ini)/abs(cost_ini);
                
                obj.optimizer_unconstr.kappa = obj.optimizer_unconstr.kappa/obj.optimizer_unconstr.kfrac;
                obj.stop_Criteria_x = ~((incr_cost < 0 && incr_vol_ls < obj.max_constr_change) || obj.optimizer_unconstr.kappa <= obj.kappa_min);
            end
            active_constr = obj.penalty > 0;
            obj.stop_criteria = obj.optimizer_unconstr.opt_cond >= obj.optimizer_unconstr.optimality_tol || any(abs(constraint.value(active_constr)) > obj.optimizer_unconstr.constr_tol(active_constr));
        end
        function checkInitial(obj)
            if isempty(obj.optimizer_unconstr.Ksmooth) 
            obj.optimizer_unconstr.Msmooth=obj.Msmooth;
            obj.optimizer_unconstr.Ksmooth=obj.Ksmooth;
            end
        end
        
    end
    
end