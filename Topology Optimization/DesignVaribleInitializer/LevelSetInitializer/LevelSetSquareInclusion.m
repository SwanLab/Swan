classdef LevelSetSquareInclusion <  ...
         LevelSetAbstractSquareInclusion & ...
         MaxNormDescriptor
        
    methods (Access = public)
        
        function obj = LevelSetSquareInclusion(cParams)
            obj.m = cParams.widthSquare;
            obj.compute(cParams);
        end
    end    
  
end

