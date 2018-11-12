classdef DesignVaribleInitializer_Feasible < DesignVaribleInitializer
    methods
        function obj = DesignVaribleInitializer_Feasible(settings,mesh,epsilon)
            obj@DesignVaribleInitializer(settings,mesh,epsilon);
        end
        
        function x = compute_initial_x(obj)
            initial_holes = false(size(obj.mesh.coord,1),1);
            obj.x(initial_holes) = obj.hole_value;
            x = obj.x;
        end
    end
end

