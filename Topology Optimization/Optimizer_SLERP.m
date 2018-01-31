classdef Optimizer_SLERP < Optimizer
    properties
        kfrac
        max_constr_change
        volume_initial
        opt_cond
        kappa_min
        optimality_tol
        constr_tol
        nconstr
    end 
    methods
        function obj=Optimizer_SLERP(settings)
            obj@Optimizer(settings);
            obj.kappa=1;
            obj.kappa_min=1e-15;
            obj.max_constr_change=+Inf;
            obj.kfrac=2;
            obj.nconstr=settings.nconstr;
        end
        function optimality_tol=get.optimality_tol(obj)
            optimality_tol=0.0175*obj.target_parameters.optimality_tol/1e-3;
        end
        function constr_tol=get.constr_tol(obj)
            constr_tol(1:obj.nconstr)=obj.target_parameters.constr_tol;
        end
        function x_ls=updateX(obj,x_ini,cost,constraint,interpolation,filter)
                x_ls=obj.updatePhi(x_ini,obj.objfunc.gradient);
                obj.update_physical_variables(x_ls,interpolation,filter);                
                cost.computef(x_ls,obj.physicalProblem,interpolation,filter);
                constraint.computef(x_ls,obj.physicalProblem,interpolation,filter);
                obj.shfunc_volume.computef(x_ls,obj.physicalProblem,interpolation,filter);
                
                obj.objfunc.computeFunction(cost,constraint)
                volume_ls =obj.shfunc_volume.value;
                
                incr_vol_ls = abs(volume_ls - obj.volume_initial);
                incr_cost = (obj.objfunc.value - obj.objfunc.value_initial)/abs(obj.objfunc.value_initial);
                
                obj.kappa = obj.kappa/obj.kfrac;
                obj.stop_criteria= ~((incr_cost < 0 && incr_vol_ls < obj.max_constr_change) || obj.kappa <= obj.kappa_min);            
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
        function phi=updatePhi(obj,design_variable,gradient)
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
