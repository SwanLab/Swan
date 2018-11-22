classdef LevelSetWidthDescriptor < handle
    
    properties (Access = protected)
        pos
        dist
    end
    
    properties (Access = protected, Abstract)
        levelSet
    end
    
    methods (Access = protected)
                
        function computeLevelSet(obj)
            obj.computeAdimensionalAndCenteredPosition()
            obj.computeDistance();
            obj.levelSet = obj.dist;
        end
        
    end
    
    methods (Access = protected, Static)
        
        function w = computeWidth(m,x)
            wx = max(x) - min(x);
            w = 0.5*m*wx;
        end
        
    end
    
    methods (Access = protected, Abstract)
        computeAdimensionalAndCenteredPosition(obj)
        computeDistance(obj)
    end
    
end

