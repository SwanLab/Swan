classdef testPerimeterRectangleQuadrilateral < testUnfittedIntegration_InternalIntegrator...
                                  & testUnfittedPerimeterRectangleIntegration
    
    properties (Access = protected)
        testName = 'test_rectangle_quadrilateral';
        analyticalValue = 8;
        meshType = 'BOUNDARY';
    end
    
end

