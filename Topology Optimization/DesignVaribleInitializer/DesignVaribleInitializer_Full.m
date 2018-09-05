classdef DesignVaribleInitializer_Full < DesignVaribleInitializer
    methods
        function obj = DesignVaribleInitializer_Full(settings,mesh,epsilon)
            obj@DesignVaribleInitializer(settings,mesh,epsilon);
        end
        
        function x = compute_initial_x(obj)
            x = obj.x;
        end
    end
end

