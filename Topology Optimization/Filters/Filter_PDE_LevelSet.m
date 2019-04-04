classdef Filter_PDE_LevelSet < Filter_PDE & Filter_LevelSet
    
    methods (Access = public)
        
        function obj = Filter_PDE_LevelSet(cParams)
           obj@ Filter_LevelSet(cParams);
        end
                
        function preProcess(obj)
            preProcess@Filter_PDE(obj)
            preProcess@Filter_LevelSet(obj)
        end
        
        function RHS = integrate_L2_function_with_shape_function(obj,x)
            F = ones(size(x));
            RHS = obj.computeRHS(x,F);
        end
        
        function RHS = integrate_function_along_facets(obj,x,F)
            RHS = obj.computeRHS(x,F);
        end
        
    end
    
end