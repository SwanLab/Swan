classdef LevelSetCircle < LevelSetSphereNdim
    
    methods (Access = public)
        
        function obj = LevelSetCircle(cParams)
            obj.fracRadius = cParams.fracRadius;
            obj.compute(cParams);
        end
    end
    
    methods (Access = protected)
        function computeLevelSetValue(obj)
            ls = obj.dist - 1;
            obj.levelSet = ls;
        end
    end
    
end
