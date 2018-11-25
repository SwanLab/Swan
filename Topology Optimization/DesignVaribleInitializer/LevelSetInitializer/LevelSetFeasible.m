classdef LevelSetFeasible < LevelSetCreator
    
    properties 
    end 
    
    methods (Access = public)
        
        function obj = LevelSetFeasible(input)
            obj.compute(input);            
        end
    end
    
    methods (Access = protected)
        
        function computeInitialLevelSet(obj)
            obj.computeLevelSet();
        end
        
    end
    
    methods (Access = private)
        
        function computeLevelSet(obj)
            x = obj.nodeCoord.x;
            obj.levelSet = ones(size(x,1),1);
        end
        
        
    end
    
end

