classdef LevelSetFull < LevelSetCreator
    
    methods (Access = public)
        
        function obj = LevelSetFull(cParams)
            obj.compute(cParams);
        end
        
    end
    
    methods (Access = protected)

        function computeLevelSet(obj)
           obj.levelSet = -1*ones(obj.lsSize); 
        end    
    end
    
end

