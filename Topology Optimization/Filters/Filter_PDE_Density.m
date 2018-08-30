classdef Filter_PDE_Density < Filter_PDE
    properties
    end
    
    methods
        function obj = Filter_PDE_Density(problemID,scale)
            obj@Filter_PDE(problemID,scale);
        end
        
        function rhs = integrate_L2_function_with_shape_function(obj,x)
            rhs = obj.diffReacProb.element.M*x;
        end
    end
end