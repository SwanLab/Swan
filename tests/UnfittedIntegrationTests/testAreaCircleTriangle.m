classdef testAreaCircleTriangle < testUnfittedIntegration_ExternalIntegrator...
                                  & testUnfittedSurfaceIntegration
    
    properties (Access = protected)
        testName = 'test_circle_triangle';
        analiticalArea = pi;
        meshType = 'INTERIOR';
    end
    
end

