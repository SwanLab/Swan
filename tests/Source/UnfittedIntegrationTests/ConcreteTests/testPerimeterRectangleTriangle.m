classdef testPerimeterRectangleTriangle < testUnfittedIntegration_ExternalIntegrator...
                                  & testUnfittedPerimeterRectangleIntegration
    
    properties (Access = protected)
        testName = 'test_rectangle_triangle';
        analyticalValue = 6;
        meshType = 'BOUNDARY';
    end
    
end

