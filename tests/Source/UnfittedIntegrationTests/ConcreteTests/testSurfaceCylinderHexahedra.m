classdef testSurfaceCylinderHexahedra < testUnfittedIntegration_ExternalIntegrator...
                                  & testUnfittedSurfaceCylinderIntegration
    
   properties (Access = protected)
        testName = 'test_cylinder_hexahedra';  
        analyticalValue = pi*2 + 2*pi*2;
        meshType = 'BOUNDARY';
   end
 
end

