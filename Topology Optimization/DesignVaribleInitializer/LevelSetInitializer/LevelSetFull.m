classdef LevelSetFull < LevelSetCreator
    
    methods (Access = public)
        
        function obj = LevelSetFull(input)
            obj.compute(input);
        end
        
    end
    
    methods (Access = protected)
        
        function computeInitialLevelSet(obj)
            obj.createLevelSet()
            obj.computeDesignVariable()
        end
        
    end
    
    methods (Access = private)
        function createLevelSet(obj)
           obj.levelSet = -1*ones(obj.lsSize); 
        end
               
    end
    
end

