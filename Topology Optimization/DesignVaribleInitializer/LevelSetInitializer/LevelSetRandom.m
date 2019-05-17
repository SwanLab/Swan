classdef LevelSetRandom < LevelSetCreator
    
    methods (Access = public)
        
        function obj = LevelSetRandom(cParams)
            obj.compute(cParams);
        end
        
    end
    
    methods (Access = protected)
        
        function computeLevelSet(obj)
            x = obj.nodeCoord;
            obj.levelSet = 2*rand(size(x,1),1) - 1;
        end              
        
    end
    
end

