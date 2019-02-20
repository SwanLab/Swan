classdef LevelSetSquareInclusion <  ...
         LevelSetAbstractSquareInclusion & ...
         MaxNormDescriptor
        
    methods (Access = public)
        
        function obj = LevelSetSquareInclusion(input)
            obj.m = input.widthSquare;
            obj.compute(input);
        end
    end    
  
end

