classdef DesignVaribleInitializer_Feasible < LevelSetCreator
    methods
        function obj = DesignVaribleInitializer_Feasible(input)
            obj.compute(input);
        end
        
    end
    
    methods (Access = protected)
        
        function x = computeInitialLevelSet(obj)
            initial_holes = false(size(obj.mesh.coord,1),1);
            obj.x(initial_holes) = obj.hole_value;
            x = obj.x;
        end
    end
end

