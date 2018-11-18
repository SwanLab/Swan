classdef DesignVaribleInitializer_Horizontal < LevelSetCreator
    properties
        
    end
    
    methods
        function obj = DesignVaribleInitializer_Horizontal(input)
            obj.compute(input);
        end
        
    end
    
    methods (Access = protected)
        
        function x = computeInitialLevelSet(obj)
%            initial_holes = obj.mesh.coord(:,2) > 0.6 | obj.mesh.coord(:,2) < 0.4;
            initial_holes = obj.mesh.coord(:,2) < 0.6 & obj.mesh.coord(:,2) > 0.4;

            obj.x(initial_holes) = obj.hole_value;
            x = obj.x;
        end
    end
end

