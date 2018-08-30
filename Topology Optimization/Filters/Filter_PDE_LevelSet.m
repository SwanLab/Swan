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
        
        function rhs = integrate_facet_with_shape_function(obj,x,F)
            rhs = obj.computeRHS_facet(x,F);
        end
        
        function S = integrate_facet_to_compute_surface(obj,x)
            S = obj.computeFacetSurface(x);
        end
    end
end