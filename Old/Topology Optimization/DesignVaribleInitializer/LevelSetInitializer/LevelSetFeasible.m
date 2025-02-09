classdef LevelSetFeasible < LevelSetCreator
    
    properties 
    end 
    
    methods (Access = public)
        
        function obj = LevelSetFeasible(cParams)
            obj.compute(cParams);            
        end
    end
    
    methods (Access = protected)
        
        function computeLevelSet(obj)
            x = obj.nodeCoord;
            obj.levelSet = ones(size(x,1),1);
        end
        
    end
    
end

