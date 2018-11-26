classdef Optimizer_SLERP < Optimizer_Unconstrained
    
    properties
        optimality_tol
        theta = 0.1
    end
    
    methods
        function obj = Optimizer_SLERP(settings,epsilon)
            obj@Optimizer_Unconstrained(settings,epsilon);
            obj.max_constr_change = +Inf;
            obj.nconstr = settings.nconstr;
        end
        
        function optimality_tol = get.optimality_tol(obj)
            optimality_tol = (0.0175/1e-3)*obj.target_parameters.optimality_tol;
        end
        
        function phi = computeX(obj,design_variable,gradient)
            phi_n = design_variable;
            norm_g = sqrt(obj.scalar_product.computeSP(gradient,gradient));            
            beta1 = sin((1-obj.line_search.kappa)*obj.theta)/sin(obj.theta);
            beta2 = sin(obj.line_search.kappa*obj.theta)/sin(obj.theta);
            phi = beta1*phi_n + beta2*gradient/norm_g;
            obj.theta = obj.computeTheta(design_variable,gradient);            
        end
        
        function theta = computeTheta(obj,phi,g)
            phi = obj.normalizeFunction(phi);
            g   = obj.normalizeFunction(g);
            phiXg = obj.scalar_product.computeSP(phi,g);
            theta = real(acos(phiXg));
            obj.opt_cond = theta;
            obj.theta = theta;
        end
        
        function x = normalizeFunction(obj,x)
            norm2 = obj.scalar_product.computeSP(x,x);
            xNorm = sqrt(norm2);
            x = x/xNorm;
        end
        
    end
end