classdef testPlotCylinderHexahedra < testPlotting_Composite
    properties (Access = protected)
        testName = 'test_cylinder_hexahedra';
        meshType = 'BOUNDARY';
        meshIncludeBoxContour = true;
    end
    
    properties (GetAccess = public, SetAccess = private)
        viewAngle = [1 1 1]
    end
end

