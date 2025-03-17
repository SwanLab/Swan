classdef LevelSetWithSphereInclusion < LevelSetSphereNdim
    
    methods (Access = public)
        
        function obj = LevelSetWithSphereInclusion(cParams)
            obj.fracRadius = cParams.fracRadius;
            obj.compute(cParams);
        end
    end
    
    methods (Access = protected)
        function computeLevelSetValue(obj)
            ls = 1 - obj.dist;
            obj.levelSet = ls;
        end
    end
    
end


