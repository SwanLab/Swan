classdef Optimizer_PG < Optimizer_Unconstrained
    properties
        kappaMultiplier
        optimality_tol
        constr_tol
    end 
    methods
        function obj = Optimizer_PG(settings,epsilon)
            obj@Optimizer_Unconstrained(settings,epsilon);
            obj.kfrac = 2;
            obj.kappaMultiplier = settings.kappaMultiplier;
            obj.kappa_min = 1e-15;
            obj.max_constr_change = +Inf;
        end
        
        function optimality_tol = get.optimality_tol(obj)
            optimality_tol = obj.target_parameters.optimality_tol;
        end
        
        function constr_tol = get.constr_tol(obj)
            constr_tol(1:obj.nconstr) = obj.target_parameters.constr_tol;
        end
        
        function x = updateX(obj,x_ini,cost,constraint)                 
                x = obj.updateRho(x_ini,obj.objfunc.gradient);
                cost.computef(x);
                constraint.computef(x);
                constraint =obj.setConstraint_case(constraint);
                
                obj.objfunc.computeFunction(cost,constraint)
%                 cost_ls = cost.value + obj.lambda*constraint.value + 0.5*obj.penalty*(constraint.value.*constraint.value);
                
                incr_norm_L2  = obj.norm_L2(x,x_ini);
                incr_cost = (obj.objfunc.value - obj.objfunc.value_initial)/abs(obj.objfunc.value_initial);
                
                obj.kappa = obj.kappa/obj.kfrac;
                obj.stop_criteria =  ~((incr_cost < 0 && incr_norm_L2 < obj.max_constr_change) || obj.kappa <=  obj.kappa_min);
                
                obj.stop_vars(1,1) = incr_cost;     obj.stop_vars(1,2) = 0;
                obj.stop_vars(2,1) = incr_norm_L2;   obj.stop_vars(2,2) = obj.max_constr_change;
                obj.stop_vars(3,1) = obj.kappa;     obj.stop_vars(3,2) = obj.kappa_min;
        end
        function rho = updateRho(obj, design_variable,gradient)  
            rho_n = design_variable;
            rho_step = rho_n-obj.kappa*gradient;
            ub = ones(length(rho_n(:,1)),1);
            lb = zeros(length(rho_n(:,1)),1);
            rho = max(min(rho_step,ub),lb);    
            obj.opt_cond = sqrt(obj.scalar_product.computeSP(rho - rho_n,rho - rho_n))/sqrt(obj.scalar_product.computeSP(rho_n,rho_n));
        end
        function computeKappa(obj,x,gradient)
            if isempty(obj.kappa)
                norm_gamma = sqrt(obj.scalar_product.computeSP(x,x));
                norm_g = sqrt(obj.scalar_product.computeSP(gradient,gradient));
                obj.kappa = norm_gamma/norm_g;
            else
                obj.kappa = obj.kappaMultiplier*obj.kappa*obj.kfrac;
            end
        end
    end
    
end
