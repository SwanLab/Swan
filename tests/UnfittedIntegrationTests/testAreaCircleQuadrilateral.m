classdef testAreaCircleQuadrilateral < testUnfittedIntegration...
                                  & testUnfittedSurfaceIntegration
    
    properties (Access = protected)
        testName = 'test_circle_quadrilateral';
        analiticalArea = pi;
        meshType = 'INTERIOR';
    end
    
end

