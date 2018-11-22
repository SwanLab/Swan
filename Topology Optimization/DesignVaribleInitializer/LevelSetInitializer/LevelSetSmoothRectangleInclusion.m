classdef LevelSetSmoothRectangleInclusion < ...
         LevelSetAbstractRectangleInclusion & ...
         PnormDescriptor
        
    methods (Access = public)
        
        function obj = LevelSetSmoothRectangleInclusion(input)
            obj.m1 = input.m1;
            obj.m2 = input.m2;
            obj.pnorm  = input.p;
            obj.compute(input);
        end
    end
    

    
end

