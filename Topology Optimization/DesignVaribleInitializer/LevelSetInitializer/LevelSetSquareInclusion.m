classdef LevelSetSquareInclusion <  ...
         LevelSetAbstractSquareInclusion & ...
         MaxNormDescriptor
        
    methods (Access = public)
        
        function obj = LevelSetSquareInclusion(input)
            obj.m = input.m;
            obj.compute(input);
        end
    end    
  
end

