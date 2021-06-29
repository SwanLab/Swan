classdef LevelSetGiven < LevelSetCreator
    
    properties (Access = private)
       value 
    end
            
    methods (Access = public)
        
        function obj = LevelSetGiven(cParams)
            obj.value = cParams.value;
            obj.compute(cParams);
        end
        
    end
    
    methods (Access = protected)

        function computeLevelSet(obj)
           obj.levelSet = obj.value; 
        end    
    end
    
end