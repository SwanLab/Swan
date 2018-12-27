classdef LevelSetWithCircleInclusion < LevelSetSphereNdim
   
    methods (Access = public)
        
        function obj = LevelSetWithCircleInclusion(input)
            obj.fracRadius = 0.4;
            obj.compute(input);
        end
    end
    
    methods (Access = protected)
        function computeLevelSetValue(obj)
            ls = 1 - obj.dist;
            obj.levelSet = ls;
        end
    end
    
   
end

