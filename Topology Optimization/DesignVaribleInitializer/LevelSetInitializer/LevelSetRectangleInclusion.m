classdef LevelSetRectangleInclusion < ...
         LevelSetAbstractRectangleInclusion & ...
         MaxNormDescriptor        

    methods (Access = public)
        function obj = LevelSetRectangleInclusion(input)
            obj.m1 = input.widthH;
            obj.m2 = input.widthV;
            obj.compute(input);
        end
    end      
end

