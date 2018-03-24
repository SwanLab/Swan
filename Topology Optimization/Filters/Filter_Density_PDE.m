classdef Filter_Density_PDE < Filter_PDE
    properties
    end
    
    methods
        function obj = Filter_Density_PDE(problemID,scale)
            obj@Filter_PDE(problemID,scale);
        end
        
        function rhs = integrate_L2_function_with_shape_function(obj,x)
            rhs = obj.Msmooth*x;
        end
    end
end