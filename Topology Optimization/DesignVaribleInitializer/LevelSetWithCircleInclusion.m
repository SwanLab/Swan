classdef LevelSetWithCircleInclusion < LevelSetCreator
    
    properties (Access = private)
        radius
        center        
    end
    
    methods (Access = public)
        
        function obj = LevelSetWithCircleInclusion(input)
            obj.compute(input);
        end
    end
    
    methods (Access = protected)
        
    function computeInitialLevelSet(obj)
        obj.computeRadius()
        obj.computeCircleCenter()
        obj.computeLevelSet()
        obj.computeDesignVariable()
    end
    
    end
    
    methods (Access = private)
        
        function computeRadius(obj)
            xLength = max(obj.nodeCoord.x) - min(obj.nodeCoord.x);
            yLength = max(obj.nodeCoord.y) - min(obj.nodeCoord.y);
            maxInteriorRadius = min(xLength/2,yLength/2);
            fracRadius = 0.4;
            obj.radius = fracRadius*maxInteriorRadius;
        end
        
        function computeCircleCenter(obj)
            x = obj.nodeCoord.x;
            y = obj.nodeCoord.y;
            obj.center(1) = 0.5*(max(x) + min(x));
            obj.center(2) = 0.5*(max(y) + min(y));            
        end
        
        function computeLevelSet(obj)
             x = obj.nodeCoord.x;
             y = obj.nodeCoord.y;
             cx = obj.center(1);
             cy = obj.center(2);             
             r = obj.radius;
             ls = ((x-cx)/r).^2 + ((y-cy)/r).^2 - 1;
             obj.levelSet = ls;
        end
        
        function computeDesignVariable(obj)
            phi = obj.levelSet;
            obj.x( phi < 0 ) = obj.hole_value;
        end
        
    end
    
end

