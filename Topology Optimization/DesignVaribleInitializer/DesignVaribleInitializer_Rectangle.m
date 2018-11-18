classdef DesignVaribleInitializer_Rectangle < LevelSetCreator
    properties
        
    end
    
    methods (Access = public)
       
        function obj = DesignVaribleInitializer_Rectangle(input)
           obj.compute(input)
        end    
        
    end
    
    methods (Access = protected)
        
        function x = computeInitialLevelSet(obj,m1,m2)
            W = max(obj.mesh.coord(:,1)) - min(obj.mesh.coord(:,1));
            H = max(obj.mesh.coord(:,2)) - min(obj.mesh.coord(:,2));
            center_x = 0.5*(max(obj.mesh.coord(:,1)) + min(obj.mesh.coord(:,1)));
            center_y = 0.5*(max(obj.mesh.coord(:,2)) + min(obj.mesh.coord(:,2)));
            offset_x = 0.5*m1*W;
            offset_y = 0.5*m2*H;
            xrange = obj.mesh.coord(:,1) < (center_x+offset_x) & obj.mesh.coord(:,1) > (center_x-offset_x);
            yrange = obj.mesh.coord(:,2) < (center_y+offset_y) & obj.mesh.coord(:,2) > (center_y-offset_y);
            initial_holes = and(xrange,yrange);
            obj.x(initial_holes) = obj.hole_value;
            x = obj.x;
        end
    end
end

