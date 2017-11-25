classdef Algorithm_PG < Algorithm
    properties
        stop_Criteria_ls=1;
        stop_Criteria_opt=1;
        lambda
        kappa
        kfrac
        max_constr_change
        kappa_min
        penalty                  
    end 
    methods
        function obj=Algorithm_PG(settings)
            obj@Algorithm(settings);
            obj.lambda=0;
            obj.optimality_tol=1e-3;
            obj.constr_tol=1e-3;
            obj.penalty=ones(settings.nconstr,1);
            obj.kfrac=4;
            obj.max_constr_change=+Inf;
            obj.kappa_min=1.0000e-15;
        end
        function updateX(obj,x_ini,cost,constraint, physProblem, interpolation,filter)  
            obj.Msmooth=physProblem.computeMass(2);
            obj.Ksmooth=physProblem.computeKsmooth;
            cost.h_C_0=cost.value;
            iter=0;
            physProblem=obj.updateEquilibrium(x_ini,physProblem,interpolation,filter);
            cost.computef(x_ini,physProblem,interpolation,filter);
            constraint.computef(x_ini,physProblem,interpolation,filter);
            obj.shfunc_volume.computef(x_ini,physProblem,interpolation,filter);
            while(obj.stop_Criteria_opt)
                iter=iter+1
                obj.plotX(x_ini,physProblem)
                volume = obj.shfunc_volume.value;                 
                obj.lambda = obj.lambda+obj.penalty*constraint.value;                
                cost_ini = cost.value + obj.lambda*constraint.value + 0.5*obj.penalty*(constraint.value.*constraint.value);
                gradient_ini = constraint.gradient*obj.lambda' + constraint.gradient*(obj.penalty'.*constraint.value) + cost.gradient;
                obj.computeKappa(x_ini,gradient_ini);
                while(obj.stop_Criteria_ls)   
                    x_ls=obj.designVariableUpdate(x_ini,obj.kappa,gradient_ini);
                    physProblem=obj.updateEquilibrium(x_ls,physProblem,interpolation,filter);
                    cost.computef(x_ls,physProblem,interpolation,filter);
                    constraint.computef(x_ls,physProblem,interpolation,filter);
                    obj.shfunc_volume.computef(x_ls,physProblem,interpolation,filter);
                    
                    cost_ls = cost.value + obj.lambda*constraint.value + 0.5*obj.penalty*(constraint.value.*constraint.value);             
                    volume_ls =obj.shfunc_volume.value; 
                    
                    incr_vol_ls = abs(volume_ls - volume);
                    incr_cost = (cost_ls - cost_ini)/abs(cost_ini);
                    
                    obj.kappa = obj.kappa/obj.kfrac;
                    obj.stop_Criteria_ls = ~((incr_cost < 0 && incr_vol_ls < obj.max_constr_change) || obj.kappa <= obj.kappa_min);              
                end
                obj.stop_Criteria_ls=1;         
                incre_x = sqrt(obj.scalar_product(x_ls - x_ini,x_ls - x_ini))/sqrt(obj.scalar_product(x_ini,x_ini))
                x_ini=x_ls;  
                active_constr = obj.penalty > 0;
                obj.stop_Criteria_opt = incre_x >= obj.optimality_tol || any(abs(constraint.value(active_constr)) > obj.constr_tol(active_constr));
            end        
        end
        function computeKappa(obj,x,gradient)
            if isempty(obj.kappa)
                norm_gamma = sqrt(obj.scalar_product(x,x));
                norm_g = sqrt(obj.scalar_product(gradient,gradient));
                obj.kappa = norm_gamma/norm_g;
            else
                obj.kappa = 10*obj.kappa*obj.kfrac;
            end
        end
        function rho=designVariableUpdate(obj,design_variable,kappa,gradient)
            rho_n = design_variable;
            rho_step = rho_n - 
            *gradient;
            ub = ones(length(rho_n(:,1)),1);
            lb = zeros(length(rho_n(:,1)),1);
            rho = max(min(rho_step,ub),lb);
            end

    end
    
end
