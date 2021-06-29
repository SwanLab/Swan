classdef LevelSetWithSphereInclusion < LevelSetSphereNdim
    
    methods (Access = public)
        
        function obj = LevelSetSphere(cParams)
            obj.fracRadius = cParams;
            obj.compute(cParams);
        end
    end
    
    methods (Access = protected)
        function computeLevelSet(obj)
            ls = 1 - obj.dist;
            obj.levelSet = ls;
        end
    end
    
end


