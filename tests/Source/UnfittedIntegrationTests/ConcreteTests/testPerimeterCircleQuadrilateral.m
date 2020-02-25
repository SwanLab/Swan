classdef testPerimeterCircleQuadrilateral < testUnfittedIntegration_ExternalIntegrator...
                                  & testUnfittedPerimeterIntegration
    
    properties (Access = protected)
        testName = 'test_circle_quadrilateral';
        analyticalValue = 2*pi;
        meshType = 'BOUNDARY';
    end
end

