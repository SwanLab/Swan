classdef Optimizer_PG < Optimizer
    properties
        kfrac
        opt_cond           
    end 
    methods
        function obj=Optimizer_PG(settings)
            obj@Optimizer(settings);
            obj.optimality_tol=1e-3;
            obj.constr_tol=1e-3;
            obj.kfrac=2;
        end
        function rho=updateX(obj, design_variable,gradient)  
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
