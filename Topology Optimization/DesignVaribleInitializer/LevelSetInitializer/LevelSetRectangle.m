classdef LevelSetRectangle < ...
        LevelSetAbstractRectangleInclusion & ...
        MaxNormDescriptor
    
    properties (Access = private)
        rectangleInclusion
    end
    
    methods (Access = public)
        
        function obj = LevelSetRectangle(cParams)
            obj.m1 = cParams.widthH;
            obj.m2 = cParams.widthV;
            obj.compute(cParams);
            obj.createRectangleInclusion(cParams);
            obj.reverseLevelSet();
        end
        
    end
    
    methods (Access = private)
        
        function createRectangleInclusion(obj,cParams)
            ls = LevelSetRectangleInclusion(cParams);
            obj.rectangleInclusion = ls;
        end
        
        function reverseLevelSet(obj)
            ls = obj.levelSet;
            obj.levelSet = -ls;
        end
        
    end
    
end

