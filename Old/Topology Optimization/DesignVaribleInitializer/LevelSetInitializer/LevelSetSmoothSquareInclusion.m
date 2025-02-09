classdef LevelSetSmoothSquareInclusion <  ...
         LevelSetAbstractSquareInclusion & ...
         PnormDescriptor
     
    methods (Access = public)
        
        function obj = LevelSetSmoothSquareInclusion(cParams)
            obj.m = cParams.widthSquare;
            obj.pnorm = cParams.pnorm;
            obj.compute(cParams);
        end
    end    
  
end

