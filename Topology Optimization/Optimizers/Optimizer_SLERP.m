classdef Optimizer_SLERP < Optimizer_Unconstrained
    
    properties
        optimality_tol
    end
    
    methods
        function obj = Optimizer_SLERP(settings,epsilon)
            obj@Optimizer_Unconstrained(settings,epsilon);
            obj.ini_design_value = -1.015243959022692;
            obj.hole_value = 0.507621979511346;
            obj.kappa = 1;
            obj.kappa_min = 1e-5;
            obj.max_constr_change = +Inf;
            obj.kfrac = 2;
            obj.nconstr = settings.nconstr;
        end
        
        function optimality_tol = get.optimality_tol(obj)
            optimality_tol = (0.0175/1e-3)*obj.target_parameters.optimality_tol;
        end
        
        function phi = computeX(obj,design_variable,gradient)
            theta = obj.computeTheta(design_variable,gradient);
            phi_n = design_variable;
            norm_g = sqrt(obj.scalar_product.computeSP(gradient,gradient));
            
            beta1 = sin((1-obj.kappa)*theta)/sin(theta);
            beta2 = sin(obj.kappa*theta)/sin(theta);
            phi = beta1*phi_n + beta2*gradient/norm_g;
        end
        
        function theta = computeTheta(obj,phi,g)
            norm_phi = sqrt(obj.scalar_product.computeSP(phi,phi));
            norm_g = sqrt(obj.scalar_product.computeSP(g,g));
            %norm_g_f = sqrt(scalar_product(g/norm_g -phi,g/norm_g -phi));
            scl_phi_g = obj.scalar_product.computeSP(phi,g);
            theta = real(acos(scl_phi_g/(norm_phi*norm_g)));
            obj.opt_cond = theta;
            %norm_dif_rel = norm_g_f;
        end
        
        function initKappa(obj,~,~,~)
            obj.kappa = 1;
        end
        
        function computeKappa(obj)
            obj.kappa = obj.kappa/obj.kfrac;
        end
    end
end
