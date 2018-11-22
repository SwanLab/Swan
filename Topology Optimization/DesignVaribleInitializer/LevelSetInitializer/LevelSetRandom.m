classdef LevelSetRandom < LevelSetCreator
    
    properties 
    end
    
    methods (Access = protected)
        
        function computeInitialLevelSet(obj)
            obj.computeLevelSet()
            %obj.computeDesignVariable()            
        end
                       
    end
    
    methods (Access = private)
        
        function computeLevelSet(obj)
            x = obj.nodeCoord.x;
            obj.levelSet = 2*rand(size(x,1),1) - 1;
        end
       
%         function computeDesignVariable(obj)
%             phi = obj.levelSet;
%             obj.x(phi < 0) = obj.hole_value;
%         end        
        
    end
    
end

