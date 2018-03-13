classdef Optimizer_Unconstrained < handle
    %Optimizer_Unconstrained Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        scalar_product
        objfunc
        target_parameters
        kfrac
        max_constr_change
        opt_cond
        kappa
        kappa_min
        %         optimality_tol
        %         constr_tol
        nconstr
        stop_vars
        stop_criteria = 1;
    end
    
    methods
        function obj = Optimizer_Unconstrained(settings,epsilon)
            obj.scalar_product = ScalarProduct(settings.filename,epsilon);
        end
    end
    
    methods (Access = protected)
        function N_L2 = norm_L2(obj,x,x_ini)
            inc_x = x-x_ini;
            N_L2 = obj.scalar_product.computeSP_M(inc_x,inc_x)/obj.scalar_product.computeSP_M(x_ini,x_ini);
        end
    end
end

