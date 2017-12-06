classdef Optimizer_PG < Optimizer
    properties
        kfrac
        opt_cond  
        volume_initial
        kappa_min
        max_constr_change
    end 
    methods
        function obj=Optimizer_PG(settings)
            obj@Optimizer(settings);
            obj.optimality_tol=1e-3;
            obj.constr_tol(1:settings.nconstr)=1e-3;
            obj.kfrac=2;
            obj.kappa_min=1e-15;
            obj.max_constr_change=+Inf;
        end
        function x=updateX(obj,x_ini,cost,constraint,physProblem,interpolation,filter) 
                x=obj.updateRho(x_ini,obj.objfunc.gradient);
                physProblem=obj.updateEquilibrium(x,physProblem,interpolation,filter);
                cost.computef(x,physProblem,interpolation,filter);
                constraint.computef(x,physProblem,interpolation,filter);
                obj.shfunc_volume.computef(x,physProblem,interpolation,filter);
                
                obj.objfunc.computeFunction(cost,constraint)
%                 cost_ls = cost.value + obj.lambda*constraint.value + 0.5*obj.penalty*(constraint.value.*constraint.value);
                volume_ls =obj.shfunc_volume.value;
                
                incr_vol_ls = abs(volume_ls - obj.volume_initial);
                incr_cost = (obj.objfunc.value - obj.objfunc.value_initial)/abs(obj.objfunc.value_initial);
                
                obj.kappa = obj.kappa/obj.kfrac;
                obj.stop_criteria= ~((incr_cost < 0 && incr_vol_ls < obj.max_constr_change) || obj.kappa <= obj.kappa_min);  
        end
        function rho=updateRho(obj, design_variable,gradient)  
            rho_n = design_variable;
            rho_step = rho_n-obj.kappa*gradient;
            ub = ones(length(rho_n(:,1)),1);
            lb = zeros(length(rho_n(:,1)),1);
            rho = max(min(rho_step,ub),lb);    
            obj.opt_cond = sqrt(obj.scalar_product(rho - rho_n,rho - rho_n))/sqrt(obj.scalar_product(rho_n,rho_n));
        end
        function computeKappa(obj,x,gradient)
            if isempty(obj.kappa)
                norm_gamma = sqrt(obj.scalar_product(x,x));
                norm_g = sqrt(obj.scalar_product(gradient,gradient));
                obj.kappa = norm_gamma/norm_g;
            else
                obj.kappa = 1*obj.kappa*obj.kfrac;
            end
        end
    end
    
end
