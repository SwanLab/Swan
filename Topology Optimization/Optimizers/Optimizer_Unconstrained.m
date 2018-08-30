classdef Optimizer_Unconstrained < Optimizer
    %Optimizer_Unconstrained Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        constr_tol
        scalar_product
        objfunc
        kfrac
        max_constr_change
        opt_cond
        good_design
        kappa
        kappa_min
    end
    
    methods
        function obj = Optimizer_Unconstrained(settings,epsilon)
            obj@Optimizer(settings);
            
            %  !! Currently only used for Unconstrained_Optimizers !!
            % (Move to Optimizer when having not black-box Constrained_Optimizer(s))
            obj.scalar_product = ScalarProduct(settings.filename,epsilon);
        end
        
        function x = updateX(obj,x_ini,cost,constraint)
            x = obj.computeX(x_ini,obj.objfunc.gradient);
            cost.computef(x);
            constraint.computef(x);
            constraint = obj.setConstraint_case(constraint);
            obj.objfunc.computeFunction(cost,constraint)
            
            incr_norm_L2  = obj.norm_L2(x,x_ini);
            incr_cost = (obj.objfunc.value - obj.objfunc.value_initial)/abs(obj.objfunc.value_initial);
            
            obj.good_design = incr_cost < 0 && incr_norm_L2 < obj.max_constr_change;
            obj.has_converged = obj.good_design || obj.kappa <= obj.kappa_min;
            
            obj.stop_vars(1,1) = incr_cost;     obj.stop_vars(1,2) = 0;
            obj.stop_vars(2,1) = incr_norm_L2;  obj.stop_vars(2,2) = obj.max_constr_change;
            obj.stop_vars(3,1) = obj.kappa;     obj.stop_vars(3,2) = obj.kappa_min;
            
            if ~obj.has_converged
                obj.computeKappa;
            end
        end
        
        function constr_tol = get.constr_tol(obj)
            constr_tol(1:obj.nconstr) = obj.target_parameters.constr_tol;
        end
    end
    
    methods (Access = public)
        function N_L2 = norm_L2(obj,x,x_ini)
            inc_x = x-x_ini;
            N_L2 = obj.scalar_product.computeSP_M(inc_x,inc_x)/obj.scalar_product.computeSP_M(x_ini,x_ini);
        end
    end
end