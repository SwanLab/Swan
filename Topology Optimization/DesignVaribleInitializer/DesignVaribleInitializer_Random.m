classdef DesignVaribleInitializer_Random < DesignVaribleInitializer
    methods
        function obj = DesignVaribleInitializer_Random(settings,mesh,epsilon)
            obj@DesignVaribleInitializer(settings,mesh,epsilon);
        end
        
        function x = compute_initial_x(obj)
            initial_holes = rand(size(obj.mesh.coord,1),1) > 0.1;
            obj.x(initial_holes) = obj.hole_value;
            x = obj.x;
        end
    end
end

