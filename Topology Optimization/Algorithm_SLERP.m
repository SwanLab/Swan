classdef Algorithm_SLERP < Algorithm
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
        function obj=Algorithm_SLERP(settings)
            obj@Algorithm(settings);
            obj.kappa=1;
            obj.lambda=0;
            obj.penalty=ones(settings.nconstr,1);
            obj.kfrac=2;
            obj.max_constr_change=+Inf;
            obj.kappa_min=1.0000e-15;
        end
        function updateX(obj,x_ini,cost,constraint, physProblem, interpolation,filter)  
            obj.Msmooth=physProblem.computeMass(2);
            obj.Ksmooth=physProblem.computeKsmooth;
            cost.h_C_0=cost.value;
            iter=0;
            cost.computef(x_ini,physProblem,interpolation,filter);
            constraint.computef(x_ini,physProblem,interpolation,filter);
            while(obj.stop_Criteria_opt)
                iter=iter+1;
                obj.plotX(x_ini,physProblem)
                volume = constraint.value;                
                obj.lambda = obj.lambda+obj.penalty*constraint.value;                
                cost_ini = cost.value + obj.lambda*constraint.value + 0.5*obj.penalty*(constraint.value.*constraint.value);
                gradient_ini = constraint.gradient*obj.lambda' + constraint.gradient*(obj.penalty'.*constraint.value) + cost.gradient;
                theta = obj.computeTheta(x_ini,gradient_ini);                
                while(obj.stop_Criteria_ls)                     
                    x_ls=obj.designVariableUpdate(x_ini,obj.kappa,theta,gradient_ini);  
                    cost.computef(x_ls,physProblem,interpolation,filter);
                    constraint.computef(x_ls,physProblem,interpolation,filter);
                    
                    cost_ls = cost.value + obj.lambda*constraint.value + 0.5*obj.penalty*(constraint.value.*constraint.value);             
                    volume_ls =constraint.value; %obj.shfunc_volume(x_ls,physProblem,interpolation,filter) 
                    gradient_ls=constraint.gradient*obj.lambda' + constraint.gradient*(obj.penalty'.*constraint.value) + cost.gradient;
                    theta_ls = obj.computeTheta(x_ls,gradient_ls);
                    incr_vol_ls = abs(volume_ls - volume);
                    incr_cost = (cost_ls - cost_ini)/abs(cost_ini);
                    
                    obj.kappa = obj.kappa/obj.kfrac;
                    obj.stop_Criteria_ls = ~((incr_cost < 0 && incr_vol_ls < obj.max_constr_change) || obj.kappa <= obj.kappa_min);              
                end
                obj.stop_Criteria_ls=1;
                x_ini=x_ls;                
                obj.kappa=1;
                active_constr = obj.penalty > 0;
                obj.stop_Criteria_opt = theta_ls >= obj.optimality_tol || any(abs(constraint.value(active_constr)) > obj.constr_tol(active_constr));
            end        
        end
        function theta=computeTheta(obj,phi,g)
            norm_phi = sqrt(obj.scalar_product(phi,phi));
            norm_g = sqrt(obj.scalar_product(g,g));
            %norm_g_f = sqrt(scalar_product(g/norm_g -phi,g/norm_g -phi));
            scl_phi_g = obj.scalar_product(phi,g);
            theta = real(acos(scl_phi_g/(norm_phi*norm_g)));
            %norm_dif_rel = norm_g_f;
        end
        function phi=designVariableUpdate(obj,design_variable,kappa,theta,gradient)
            phi_n = design_variable;
            norm_g = sqrt(obj.scalar_product(gradient,gradient));
            
            beta1 = sin((1-kappa)*theta)/sin(theta);
            beta2 = sin(kappa*theta)/sin(theta);
            phi = beta1*phi_n + beta2*gradient/norm_g;
        end

    end
    
end
