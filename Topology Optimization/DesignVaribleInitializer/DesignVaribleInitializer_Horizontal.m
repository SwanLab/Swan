classdef DesignVaribleInitializer_Horizontal < DesignVaribleInitializer
    properties
        
    end
    
    methods
        function obj = DesignVaribleInitializer_Horizontal(settings,mesh,epsilon)
            obj@DesignVaribleInitializer(settings,mesh,epsilon);
        end
        
        function x = compute_initial_x(obj)
            initial_holes = obj.mesh.coord(:,2) > 0.6 | obj.mesh.coord(:,2) < 0.4;
            obj.x(initial_holes) = obj.hole_value;
            x = obj.x;
        end
    end
end

