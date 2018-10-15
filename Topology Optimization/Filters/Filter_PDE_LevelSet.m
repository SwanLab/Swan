classdef Filter_PDE_LevelSet < Filter_PDE
    methods (Abstract)
        preProcess(obj)
    end
    
    methods
        function obj = Filter_PDE_LevelSet(problemID,scale)
            obj@Filter_PDE(problemID,scale);
        end
        
        function rhs = integrate_L2_function_with_shape_function(obj,x)
            rhs = obj.computeRHS(x);
        end
        
        function rhs = integrate_function_along_facets(obj,x,F)
            rhs = obj.computeRHS(x,F);
        end
    end
end