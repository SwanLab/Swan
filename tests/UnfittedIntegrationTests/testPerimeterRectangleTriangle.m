classdef testPerimeterRectangleTriangle < testUnfittedIntegration...
                                  & testUnfittedPerimeterRectangleIntegration
    
    properties (Access = protected)
        testName = 'test_rectangle_triangle';
        analiticalArea = 8;
        meshType = 'BOUNDARY';
    end
end

