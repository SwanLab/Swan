classdef LevelSetRectangleInclusion < ...
         LevelSetAbstractRectangleInclusion & ...
         MaxNormDescriptor        

    methods (Access = public)
        
        function obj = LevelSetRectangleInclusion(input)
            obj.m1 = input.m1;
            obj.m2 = input.m2;
            obj.compute(input);
        end
    end
          
    
end

