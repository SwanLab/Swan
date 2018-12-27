classdef Filter_PDE_Density < Filter_PDE
    methods (Access = public)        
        function RHS = integrate_L2_function_with_shape_function(obj,x)
            RHS = obj.diffReacProb.element.M*x;
        end
    end
end