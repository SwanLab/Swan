classdef LevelSetWithCircleInclusion < LevelSetSphereNdim
   
    methods (Access = public)
        
        function obj = LevelSetWithCircleInclusion(input)
            obj.fracRadius = 0.4;
            obj.compute(input);
        end
    end
    
   
end

