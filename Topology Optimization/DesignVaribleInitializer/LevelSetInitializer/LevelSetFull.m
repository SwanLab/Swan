classdef LevelSetFull < LevelSetCreator
    
    methods (Access = public)
        
        function obj = LevelSetFull(input)
            obj.compute(input);
        end
        
    end
    
    methods (Access = protected)

        function computeLevelSet(obj)
           obj.levelSet = -1*ones(obj.lsSize); 
        end    
    end
    
end

