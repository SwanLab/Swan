classdef testPlotSphereTetrahedra < testPlotting
    properties (Access = protected)
        testName = 'test_sphere_tetrahedra';
        meshType = 'BOUNDARY';
        meshIncludeBoxContour = false;
    end
    
    properties (GetAccess = public, SetAccess = private)
        viewAngle = [1 1 1]
    end
end

