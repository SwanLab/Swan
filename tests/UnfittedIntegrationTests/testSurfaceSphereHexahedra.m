classdef testSurfaceSphereHexahedra < testUnfittedIntegration_InternalIntegrator...
                                  & testUnfittedSurfaceIntegration
    
    properties (Access = protected)
        testName = 'test_sphere_hexahedra';
        analiticalArea = 4*pi;
        meshType = 'BOUNDARY';
    end
    
end

