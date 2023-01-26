classdef LevelSetCylinderBuckling < LevelSetCreator
    
    properties (Access = private)
        x0
        y0
    end

    properties (Access = private)
        Radius
    end
    
    methods (Access = public)
        
        function obj = LevelSetCylinderBuckling(cParams)
            obj.Radius = cParams.Radius;
            obj.compute(cParams);
        end
    end
    
    methods (Access = protected)
        function computeLevelSet(obj)
            R = obj.Radius;
            x = obj.nodeCoord(:,1);
            y = obj.nodeCoord(:,2);
            obj.computeCylinderCentre(x,y);
            ls = (x-obj.x0).^2 + (y-obj.y0).^2 - R.^2;
            obj.levelSet = ls;
        end

        function computeCylinderCentre(obj,x,y)
            obj.x0 = 0.5*(max(x)+min(x));
            obj.y0 = 0.5*(max(y)+min(y));
        end
    end
    
end