classdef Filter_PDE_LevelSet < Filter_PDE & Filter_LevelSet
    methods
        function obj = Filter_PDE_LevelSet(problemID,scale)
            obj@Filter_PDE(problemID,scale);
        end
        
        function preProcess(obj)
            preProcess@Filter_PDE(obj)
            preProcess@Filter_LevelSet(obj)
        end
        
        function rhs = integrate_L2_function_with_shape_function(obj,x)
            rhs = obj.computeRHS(x);
        end
        
        function rhs = integrate_contour_with_shape_function(obj,x,V)
            rhs = obj.computeRHS_contour(x,V);
        end
    end
end