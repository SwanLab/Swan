classdef testAreaCircleTriangle < testUnfittedIntegration_ExternalIntegrator...
                                  & testUnfittedSurfaceIntegration
    
    properties (Access = protected)
        testName = 'test_circle_triangle';
        analyticalValue = pi;
        meshType = 'INTERIOR';
    end
    
end

