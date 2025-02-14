classdef LevelSetCylinder < LevelSetCylinderDistance
    
    methods (Access = public)
        
        function obj = LevelSetCylinder(cParams)
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


