classdef Optimizer_Unconstrained < Optimizer
    %Optimizer_Unconstrained Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        scalar_product
        objfunc
        kfrac
        max_constr_change
        opt_cond
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
    end
    
    methods (Access = protected)
        function N_L2 = norm_L2(obj,x,x_ini)
            inc_x = x-x_ini;
            N_L2 = obj.scalar_product.computeSP_M(inc_x,inc_x)/obj.scalar_product.computeSP_M(x_ini,x_ini);
        end
    end
end

