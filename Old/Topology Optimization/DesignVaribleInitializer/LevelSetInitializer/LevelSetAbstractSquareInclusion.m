classdef LevelSetAbstractSquareInclusion < ...
         LevelSetCreator & ...
         LevelSetCenterDescriptor & ...
         LevelSetWidthDescriptor
         
         
     properties (Access = protected)     
        width
        m
        pos
     end
     
     properties (Access = protected, Abstract)
        dist 
     end
    
    methods (Access = protected)
        
        function computeLevelSet(obj)
            obj.computeCenter();
            obj.computeInclusionWidth();
            obj.computeAdimensionalAndCenteredPosition()
            obj.computeDistance();            
            obj.computeLevelSetValue();           
        end        
        
        function computeInclusionWidth(obj)
            x = obj.nodeCoord(:,1);
            y = obj.nodeCoord(:,2);
            widthX = obj.computeWidth(obj.m,x);
            widthY = obj.computeWidth(obj.m,y);
            isEqualWidth = abs(widthX - widthY) < 1e-14;
            assert(isEqualWidth);
            obj.width = 0.5*(widthX + widthY);
        end
                
        function computeAdimensionalAndCenteredPosition(obj)
            x0 = obj.nodeCoord(:,1);
            y0 = obj.nodeCoord(:,2);
            x = x0 - obj.center(1);
            y = y0 - obj.center(2);
            w = obj.width;
            obj.pos = [x/w,y/w];
        end
        
        function computeLevelSetValue(obj)
            obj.levelSet = 1 - (obj.dist + 1e-14);
        end
               
    end
    
    methods (Access = protected, Abstract)
        computeDistance(obj)
    end
    
end

