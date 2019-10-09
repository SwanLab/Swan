classdef LevelSetRectangleInclusion < ...
         LevelSetAbstractRectangleInclusion & ...
         MaxNormDescriptor        

    methods (Access = public)
        
        function obj = LevelSetRectangleInclusion(cParams)
            obj.m1 = cParams.widthH;
            obj.m2 = cParams.widthV;
            obj.compute(cParams);
        end
        
    end    
    
end

