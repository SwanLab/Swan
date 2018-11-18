classdef DesignVaribleInitializerRoundedRectangle < LevelSetCreator
    properties
        
    end
    
    methods (Access= public)
       
        function obj = DesignVaribleInitializerRoundedRectangle()            
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
            
            x = obj.mesh.coord(:,1) - center_x;
            y = obj.mesh.coord(:,2) - center_y;
            
            p = 4;% + 20*obj.mesh.mean_cell_size;
            rp = ((x/offset_x).^p + (y/offset_y).^p).^(1/p) - 1;

            initial_holes = rp > 0;
            obj.x(initial_holes) = obj.hole_value;
            x = obj.x;
        end
    end
end

