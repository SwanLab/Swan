classdef testVolumeSphereHexahedra < testUnfittedIntegration_InternalIntegrator...
                                  & testUnfittedVolumeIntegration
    
   properties (Access = protected)
        testName = 'test_sphere_hexahedra';  
        analyticalValue = (4/3)*pi;
        meshType = 'INTERIOR';
   end
 
end

