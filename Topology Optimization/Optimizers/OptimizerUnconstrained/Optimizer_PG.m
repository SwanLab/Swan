classdef Optimizer_PG < Optimizer_Unconstrained
    
    properties  (GetAccess = public, SetAccess = private)
        optimality_tol
    end
    
    methods (Access = public)
        
        function obj = Optimizer_PG(settings)
            obj@Optimizer_Unconstrained(settings);
            obj.maxIncrNormX = +Inf;
        end
        
        function rho = compute(obj,design_variable,gradient)
            rho_n = design_variable;
            rho_step = rho_n-obj.line_search.kappa*gradient;
            ub = ones(length(rho_n(:,1)),1);
            lb = zeros(length(rho_n(:,1)),1);
            rho = max(min(rho_step,ub),lb);
            obj.opt_cond = sqrt(obj.scalar_product.computeSP(rho - rho_n,rho - rho_n))/sqrt(obj.scalar_product.computeSP(rho_n,rho_n));
        end
        
    end
    
    methods
        
        function optimality_tol = get.optimality_tol(obj)
            optimality_tol = obj.target_parameters.optimality_tol;
        end
        
    end
    
end