classdef Optimizer_Unconstrained < Optimizer
    %Optimizer_Unconstrained Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        constr_tol
        scalar_product
        objfunc
        max_constr_change
        opt_cond
        good_design
        line_search
    end
    
    methods
        function obj = Optimizer_Unconstrained(settings,epsilon)
            obj@Optimizer(settings);
            obj.line_search = LineSearch.create(settings,epsilon);
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
            obj.stop_updating = obj.good_design || obj.line_search.kappa <= obj.line_search.kappa_min;
            
            obj.stop_vars(1,1) = incr_cost;                 obj.stop_vars(1,2) = 0;
            obj.stop_vars(2,1) = incr_norm_L2;              obj.stop_vars(2,2) = obj.max_constr_change;
            obj.stop_vars(3,1) = obj.line_search.kappa;     obj.stop_vars(3,2) = obj.line_search.kappa_min;
            
            if ~obj.stop_updating
                obj.line_search.computeKappa;
            end
        end
        
        function constr_tol = get.constr_tol(obj)
            constr_tol(1:obj.nconstr) = obj.target_parameters.constr_tol;
        end
    end
    
    methods (Access = protected)
        function N_L2 = norm_L2(obj,x,x_ini)
            inc_x = x-x_ini;
            N_L2 = obj.scalar_product.computeSP_M(inc_x,inc_x)/obj.scalar_product.computeSP_M(x_ini,x_ini);
        end
    end
end

