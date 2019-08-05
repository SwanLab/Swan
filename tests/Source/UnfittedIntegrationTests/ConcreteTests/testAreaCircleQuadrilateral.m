classdef testAreaCircleQuadrilateral < testUnfittedIntegration_ExternalIntegrator...
                                  & testUnfittedSurfaceIntegration
    
    properties (Access = protected)
        testName = 'test_circle_quadrilateral';
        analyticalValue = pi;
        meshType = 'INTERIOR';
    end
    
end

