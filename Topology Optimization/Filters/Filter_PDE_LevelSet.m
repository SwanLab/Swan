classdef Filter_PDE_LevelSet < Filter_PDE
    methods (Abstract)
        preProcess(obj)
    end
    
    methods
        function rhs = integrate_L2_function_with_shape_function(obj,x)
            obj.unfitted_mesh.computeMesh(x);
            rhs = obj.computeRHS;
        end
        
        function rhs = integrate_function_along_facets(obj,x,F)
            obj.unfitted_mesh.computeMesh(x);
            obj.unfitted_mesh.computeGlobalConnectivities;
            rhs = obj.computeRHS(F);
        end
    end
end