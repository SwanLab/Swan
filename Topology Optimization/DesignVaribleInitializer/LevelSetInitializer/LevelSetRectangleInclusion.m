classdef LevelSetRectangleInclusion < ...
         LevelSetAbstractRectangleInclusion & ...
         MaxNormDescriptor        

    methods (Access = public)
        function obj = LevelSetRectangleInclusion(cParams)
            obj.m1 = cParams.geomParams.widthH;
            obj.m2 = cParams.geomParams.widthV;
            obj.compute(cParams);
        end
    end      
end

