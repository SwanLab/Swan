classdef LevelSetWithSphereInclusion < LevelSetSphereNdim
    
    methods (Access = public)
        
        function obj = LevelSetSphere(input)
            obj.fracRadius = 1-1e-6;
            obj.compute(input);
        end
    end
    
    methods (Access = protected)
        function computeLevelSet(obj)
            ls = 1 - obj.dist;
            obj.levelSet = ls;
        end
    end
    
end


