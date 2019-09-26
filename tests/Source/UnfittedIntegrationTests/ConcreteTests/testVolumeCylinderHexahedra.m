classdef testVolumeCylinderHexahedra < testUnfittedIntegration_ExternalIntegrator...
                                  & testUnfittedVolumeIntegration
    
   properties (Access = protected)
        testName = 'test_cylinder_hexahedra';  
        analyticalValue = pi*2;
        meshType = 'INTERIOR';
   end
 
end

