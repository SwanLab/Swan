classdef LevelSetCircle2 < LevelSetSphereNdim
    
    methods (Access = public)
        
        function obj = LevelSetCircle2(input)
            obj.fracRadius = 1.2;
            obj.compute(input);
        end
    end
    
    methods (Access = protected)
        function computeLevelSetValue(obj)
            ls = obj.dist - 1;
            obj.levelSet = ls;
        end
    end
    
end


