classdef LevelSetAbstractSquareInclusion < ...
         LevelSetCreator & ...
         LevelSetCenterDescriptor & ...
         LevelSetWidthDescriptor
         
         
     properties (Access = protected)     
        width
        m
     end
    
    methods (Access = protected)
        
        function computeInitialLevelSet(obj)
            obj.computeCenter();
            obj.computeInclusionWidth();
            obj.computeLevelSet();
            obj.computeDesignVariable();
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
               
    end
    
    

    
end

