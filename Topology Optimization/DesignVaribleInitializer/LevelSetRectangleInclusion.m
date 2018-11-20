classdef LevelSetRectangleInclusion < LevelSetCreator
    
     properties (Access = private)
        center     
        widthX
        widthY
        m1
        m2
    end
    
    methods (Access = public)
        
        function obj = LevelSetRectangleInclusion(input)
            obj.m1 = input.m1;
            obj.m2 = input.m2;
            obj.compute(input);
        end
    end
    
    methods (Access = protected)
        
        function computeInitialLevelSet(obj)
            obj.computeCenter();
            obj.computeInclusionWidths();
            obj.computeLevelSet();
            obj.computeDesignVariable();
        end
        
        function computeCenter(obj)
            x = obj.nodeCoord.x;
            y = obj.nodeCoord.y;
            obj.center(1) = 0.5*(max(x) + min(x));
            obj.center(2) = 0.5*(max(y) + min(y));
        end
        
        function computeInclusionWidths(obj)
            x = obj.nodeCoord.x;
            y = obj.nodeCoord.y;
            obj.widthX = obj.computeWidth(obj.m1,x);
            obj.widthY = obj.computeWidth(obj.m2,y);
        end
        
        function computeLevelSet(obj)
            x0 = obj.nodeCoord.x;
            y0 = obj.nodeCoord.y;
            x = x0 - obj.center(1);
            y = y0 - obj.center(2);
            wx = obj.widthX;
            wy = obj.widthY;
            ls = max(abs([x/wx,y/wy]),[],2) + 1e-14 - 1;
            obj.levelSet = ls;
        end
        
        function computeDesignVariable(obj)
            phi = obj.levelSet;
            obj.x( phi < 0) = obj.hole_value;
        end
    end
    
end

