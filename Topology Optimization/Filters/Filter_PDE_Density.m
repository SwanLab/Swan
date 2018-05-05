classdef Filter_PDE_Density < Filter_PDE
    properties
    end
    
    methods
        function rhs = integrate_L2_function_with_shape_function(obj,x)
            rhs = obj.Msmooth*x;
        end
    end
end