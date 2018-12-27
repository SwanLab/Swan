classdef LevelSetRandom < LevelSetCreator
    
    properties 
    end
    
    methods (Access = protected)
        
        function computeInitialLevelSet(obj)
            obj.computeLevelSet()
        end
                       
    end
    
    methods (Access = private)
        
        function computeLevelSet(obj)
            x = obj.nodeCoord.x;
            obj.levelSet = 2*rand(size(x,1),1) - 1;
        end              
        
    end
    
end

