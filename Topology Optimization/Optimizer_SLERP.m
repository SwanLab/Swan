classdef Optimizer_SLERP < Optimizer
    properties
        kfrac
        opt_cond
    end 
    methods
        function obj=Optimizer_SLERP(settings)
            obj@Optimizer(settings);
            obj.kappa=1;
            obj.optimality_tol=0.0175;
            obj.constr_tol=1e-3;
            obj.kfrac=2;
        end
        function theta=computeTheta(obj,phi,g)
            norm_phi = sqrt(obj.scalar_product(phi,phi));
            norm_g = sqrt(obj.scalar_product(g,g));
            %norm_g_f = sqrt(scalar_product(g/norm_g -phi,g/norm_g -phi));
            scl_phi_g = obj.scalar_product(phi,g);
            theta = real(acos(scl_phi_g/(norm_phi*norm_g)));
            obj.opt_cond=theta;
            %norm_dif_rel = norm_g_f;
        end
        function phi=updateX(obj,design_variable,gradient)
            theta=obj.computeTheta(design_variable,gradient);
            phi_n = design_variable;
            norm_g = sqrt(obj.scalar_product(gradient,gradient));
            
            beta1 = sin((1-obj.kappa)*theta)/sin(theta);
            beta2 = sin(obj.kappa*theta)/sin(theta);
            phi = beta1*phi_n + beta2*gradient/norm_g;
        end

        function computeKappa(obj,~,~,~)
            obj.kappa=1;
        end
    end
    
end
