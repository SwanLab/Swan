classdef LevelSetSphere < LevelSetSphereNdim
    
    methods (Access = public)
        
        function obj = LevelSetSphere(input)
            obj.fracRadius = 1-1e-6;
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


