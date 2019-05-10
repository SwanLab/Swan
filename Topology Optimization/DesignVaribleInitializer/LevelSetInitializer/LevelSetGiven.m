classdef LevelSetGiven < LevelSetCreator
    
    properties (Access = private)
       value 
    end
            
    methods (Access = public)
        
        function obj = LevelSetGiven(input)
            obj.value = input.value;
            obj.compute(input);
        end
        
    end
    
    methods (Access = protected)

        function computeLevelSet(obj)
           obj.levelSet = obj.value; 
        end    
    end
    
end