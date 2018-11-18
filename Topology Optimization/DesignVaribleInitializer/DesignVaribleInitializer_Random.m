classdef DesignVaribleInitializer_Random < DesignVaribleInitializer
   
    
    methods (Access = public)
        
        function obj = DesignVaribleInitializer_Random(input)
           obj.compute(input) 
        end
        
    end
    
    methods (Access = protected)
        
        function x = computeInitialLevelSet(obj)
            initial_holes = rand(size(obj.mesh.coord,1),1) > 0.1;
            obj.x(initial_holes) = obj.hole_value;
            x = obj.x;
        end
    end
end

