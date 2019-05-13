classdef LevelSetHorizontalInclusion < LevelSetCreator
    
    properties (Access = private)
    end
    
    methods (Access = public)
        
        function obj = LevelSetHorizontalInclusion(input)
            obj.compute(input);
        end
    end
    
    methods (Access = protected)
        
        function computeInitialLevelSet(obj)
            obj.computeLevelSet()
        end
        
    end
    
    methods (Access = private)
        
        function computeLevelSet(obj)
            m = 0.2;
            y0 = obj.nodeCoord.y;
            H = max(y0) - min(y0);
            center_y = 0.5*(max(y0) + min(y0));
            offset_y = 0.5*m*H;            
            y = y0 - center_y;
            obj.levelSet =  - (max(y/offset_y) - 1);
        end
        
    end
    
end

