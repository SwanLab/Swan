classdef testPerimeterRectangleTriangle < testUnfittedIntegration_ExternalIntegrator...
                                  & testUnfittedPerimeterRectangleIntegration
    
    properties (Access = protected)
        testName = 'test_rectangle_triangle';
        analyticalValue = 8;
        meshType = 'BOUNDARY';
    end
    
end

