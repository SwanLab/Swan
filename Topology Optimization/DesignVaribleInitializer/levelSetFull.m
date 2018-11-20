classdef LevelSetFull < LevelSetCreator
    
    methods (Access = public)
        
        function obj = LevelSetFull(input)
            obj.compute(input);
        end
        
    end
    
    methods (Access = protected)
        
        function computeInitialLevelSet(obj)
            obj.createLevelSet()
            obj.createDesignVariable()
        end
        
    end
    
    methods (Access = private)
        function createLevelSet(obj)
           obj.levelSet = -1*ones(obj.lsSize); 
        end
        
        function createDesignVariable(obj)
            phi = obj.levelSet;
            obj.x = obj.ini_design_value;
            obj.x(phi<0) = obj.hole_value;
        end
        
    end
    
end

