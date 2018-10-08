classdef DesignVaribleInitializer_Circle < DesignVaribleInitializer
    properties
        radius = 1-1e-6;
    end
    
    methods
        function obj = DesignVaribleInitializer_Circle(settings,mesh,epsilon)
            obj@DesignVaribleInitializer(settings,mesh,epsilon);
        end
        
        function x = compute_initial_x(obj)
            dim(1) = max(obj.mesh.coord(:,1)) - min(obj.mesh.coord(:,1));
            dim(2) = max(obj.mesh.coord(:,2)) - min(obj.mesh.coord(:,2));
            center(1) = 0.5*(max(obj.mesh.coord(:,1)) + min(obj.mesh.coord(:,1)));
            center(2) = 0.5*(max(obj.mesh.coord(:,2)) + min(obj.mesh.coord(:,2)));
            radius = obj.radius*min(dim)/2;
            phi = (obj.mesh.coord(:,1)-center(1)).^2 + (obj.mesh.coord(:,2)-center(2)).^2 - radius^2;
            
            switch obj.optimizer
                case {'SLERP','HAMILTON-JACOBI'}
                    obj.x = phi;
                otherwise
                    initial_holes = ceil(max(phi,0))>0;
                    obj.x(initial_holes) = obj.hole_value;
            end
            
            x = obj.x;
        end
    end
end

