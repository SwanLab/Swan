classdef LevelSetSmoothRectangleInclusion < ...
        LevelSetAbstractRectangleInclusion & ...
        PnormDescriptor
    
    methods (Access = public)
        
        function obj = LevelSetSmoothRectangleInclusion(cParams)
            obj.m1 = cParams.widthH;
            obj.m2 = cParams.widthV;
            obj.pnorm  = cParams.pnorm;
            obj.compute(cParams);
        end
        
    end
    
end

