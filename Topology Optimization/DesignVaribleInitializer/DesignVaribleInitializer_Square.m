classdef DesignVaribleInitializer_Square < DesignVaribleInitializer
    properties
        
    end
    
    methods
        function obj = DesignVaribleInitializer_Square(settings,mesh,epsilon)
            obj@DesignVaribleInitializer(settings,mesh,epsilon);
        end
        
        function x = compute_initial_x(obj)
            W = max(obj.mesh.coord(:,1)) - min(obj.mesh.coord(:,1));
            H = max(obj.mesh.coord(:,2)) - min(obj.mesh.coord(:,2));
            center_x = 0.5*(max(obj.mesh.coord(:,1)) + min(obj.mesh.coord(:,1)));
            center_y = 0.5*(max(obj.mesh.coord(:,2)) + min(obj.mesh.coord(:,2)));
            offset_x = 0.2*W;
            offset_y = 0.2*H;
            xrange = obj.mesh.coord(:,1) < (center_x+offset_x) & obj.mesh.coord(:,1) > (center_x-offset_x);
            yrange = obj.mesh.coord(:,2) < (center_y+offset_y) & obj.mesh.coord(:,2) > (center_y-offset_y);
            initial_holes = and(xrange,yrange);
            obj.x(initial_holes) = obj.hole_value;
            x = obj.x;
        end
    end
end

