classdef testPerimeterRectangleQuadrilateral < testUnfittedIntegration_InternalIntegrator...
                                  & testUnfittedPerimeterRectangleIntegration
    
    properties (Access = protected)
        testName = 'test_rectangle_quadrilateral';
        analyticalValue = 6;
        meshType = 'BOUNDARY';
    end
    
end

