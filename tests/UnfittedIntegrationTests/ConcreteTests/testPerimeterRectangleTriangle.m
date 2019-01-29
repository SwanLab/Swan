classdef testPerimeterRectangleTriangle < testUnfittedIntegration_ExternalIntegrator_Composite...
                                  & testUnfittedPerimeterRectangleIntegration
    
    properties (Access = protected)
        testName = 'test_rectangle_triangle';
        analyticalValue = 8;
        meshType = 'BOUNDARY';
    end
    
end

