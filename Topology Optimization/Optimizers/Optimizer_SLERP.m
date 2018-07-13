classdef Optimizer_SLERP < Optimizer_Unconstrained
    
    properties
        optimality_tol
        constr_tol
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
        
        function constr_tol = get.constr_tol(obj)
            constr_tol(1:obj.nconstr) = obj.target_parameters.constr_tol;
        end
        
        function x = updateX(obj,x_ini,cost,constraint)
            x = obj.updatePhi(x_ini,obj.objfunc.gradient);
            cost.computef(x);
            constraint.computef(x);
            constraint =obj.setConstraint_case(constraint); 
            obj.objfunc.computeFunction(cost,constraint)
            
            incr_norm_L2  = obj.norm_L2(x,x_ini);
            incr_cost = (obj.objfunc.value - obj.objfunc.value_initial)/abs(obj.objfunc.value_initial);
            
            obj.kappa = obj.kappa/obj.kfrac;
            obj.stop_criteria = ~((incr_cost < 0 && incr_norm_L2 < obj.max_constr_change) || obj.kappa <= obj.kappa_min);
            
            obj.stop_vars(1,1) = incr_cost;     obj.stop_vars(1,2) = 0;
            obj.stop_vars(2,1) = incr_norm_L2;   obj.stop_vars(2,2) = obj.max_constr_change;
            obj.stop_vars(3,1) = obj.kappa;     obj.stop_vars(3,2) = obj.kappa_min;
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
        
        function phi = updatePhi(obj,design_variable,gradient)
            theta = obj.computeTheta(design_variable,gradient);
            phi_n = design_variable;
            norm_g = sqrt(obj.scalar_product.computeSP(gradient,gradient));
            
            beta1 = sin((1-obj.kappa)*theta)/sin(theta);
            beta2 = sin(obj.kappa*theta)/sin(theta);
            phi = beta1*phi_n + beta2*gradient/norm_g;
        end
        
        function computeKappa(obj,~,~,~)
            obj.kappa = 1;
        end
    end
end
