classdef DesignVaribleInitializer_Circle < DesignVaribleInitializer
    properties
        radius = 0.4;
    end
    
    methods
        function obj = DesignVaribleInitializer_Circle(input)
            obj.compute(input);
        end
    end
    methods (Access = protected)
        
        function x = computeInitialLevelSet(obj)
            dim(1) = max(obj.mesh.coord(:,1)) - min(obj.mesh.coord(:,1));
            dim(2) = max(obj.mesh.coord(:,2)) - min(obj.mesh.coord(:,2));
            center(1) = 0.5*(max(obj.mesh.coord(:,1)) + min(obj.mesh.coord(:,1)));
            center(2) = 0.5*(max(obj.mesh.coord(:,2)) + min(obj.mesh.coord(:,2)));
            obj.radius = obj.radius*min(dim)/2;
            initial_holes = (obj.mesh.coord(:,1)-center(1)).^2 + (obj.mesh.coord(:,2)-center(2)).^2 - obj.radius^2 < 0;
            obj.x(initial_holes) = obj.hole_value;
            x = obj.x;
        end
    end
end

