classdef Filter_PDE_LevelSet < Filter_PDE
    methods (Abstract)
        preProcess(obj)
    end
    
    methods (Access = public)
        function RHS = integrate_L2_function_with_shape_function(obj,x)
            obj.unfitted_mesh.computeMesh(x);
            F = ones(size(x));
            RHS = obj.computeRHS(F);
        end
        
        function RHS = integrate_function_along_facets(obj,x,F)
            obj.unfitted_mesh.computeMesh(x);
            RHS = obj.computeRHS(F);
        end
    end
end