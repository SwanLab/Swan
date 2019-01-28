classdef testPerimeterCircleTriangle < testUnfittedIntegration_ExternalIntegrator...
                                  & testUnfittedPerimeterIntegration
    
    properties (Access = protected)
        testName = 'test_circle_triangle';
        analyticalValue = 2*pi;
        meshType = 'BOUNDARY';
    end
    
end

