classdef Filter_PDE_LevelSet < Filter_PDE
    methods (Abstract)
        preProcess(obj)
    end
    
    methods (Access = public)
        function RHS = integrate_L2_function_with_shape_function(obj,x)
            F = ones(size(x));
            RHS = obj.computeRHS(x,F);
        end
        
        function RHS = integrate_function_along_facets(obj,x,F)
            RHS = obj.computeRHS(x,F);
        end
    end
end