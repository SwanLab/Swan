classdef LevelSetSmoothSquareInclusion <  ...
         LevelSetAbstractSquareInclusion & ...
         PnormDescriptor
     
    methods (Access = public)
        
        function obj = LevelSetSmoothSquareInclusion(input)
            obj.m = input.m;
            obj.pnorm = input.p;
            obj.compute(input);
        end
    end    
  
end

