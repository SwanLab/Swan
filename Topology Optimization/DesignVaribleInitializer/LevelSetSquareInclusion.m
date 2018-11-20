classdef LevelSetSquareInclusion < LevelSetCreator
    
    properties (Access = private)
        widthX
        widthY
        widthFrac
        center
    end
    
    methods (Access = public)
        
        function obj = LevelSetSquareInclusion(input)
            input.m = 0.4;
            obj.widthFrac = input.m;
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
        
    end
    methods (Access = private)
        
        function computeCenter(obj)
            x = obj.nodeCoord(:,1);
            y = obj.nodeCoord(:,2);
            obj.center(1) = 0.5*(max(x) + min(x));
            obj.center(2) = 0.5*(max(y) + min(y));
        end
        
        function computeInclusionWidths(obj)
            x = obj.nodeCoord(:,1);
            y = obj.nodeCoord(:,2);
            obj.widthX = obj.computeWidth(obj.widthFrac,x);
            obj.widthY = obj.computeWidth(obj.widthFrac,y);
        end
        
        function computeLevelSet(obj)
            x0 = obj.nodeCoord(:,1);
            y0 = obj.nodeCoord(:,2);
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
    
    methods (Access = private, Static)
        
        function w = computeWidth(m,x)
            wx = max(x) - min(x);
            w = 0.5*m*wx;
        end
        
    end
    
    
end

