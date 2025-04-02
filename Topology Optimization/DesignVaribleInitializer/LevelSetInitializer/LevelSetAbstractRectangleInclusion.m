classdef LevelSetAbstractRectangleInclusion < ...
        LevelSetCreator & ...
        LevelSetCenterDescriptor & ...
        LevelSetWidthDescriptor
    
    properties (Access = protected)
        widthX
        widthY
        m1
        m2
        pos
    end
    
    properties (Access = protected, Abstract)
        dist
    end
    
    methods (Access = protected)
        
        function computeLevelSet(obj)
            obj.computeCenter();
            obj.computeInclusionWidths();
            obj.computeAdimensionalAndCenteredPosition()
            obj.computeDistance()
            obj.computeLevelSetValue();
        end
        
        function computeInclusionWidths(obj)
            x = obj.nodeCoord(:,1);
            y = obj.nodeCoord(:,2);
            obj.widthX = obj.computeWidth(obj.m1,x);
            obj.widthY = obj.computeWidth(obj.m2,y);
        end
        
        function computeAdimensionalAndCenteredPosition(obj)
            x0 = obj.nodeCoord(:,1);
            y0 = obj.nodeCoord(:,2);
            x = x0 - obj.center(1);
            y = y0 - obj.center(2);
            wx = obj.widthX;
            wy = obj.widthY;
            obj.pos = [x./wx,y./wy];
        end
        
        function computeLevelSetValue(obj)
            obj.levelSet = 1 - (obj.dist + 1e-14);
        end
        
    end
    
    methods (Access = protected, Abstract)
        computeDistance(obj)
    end  
end

