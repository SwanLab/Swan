classdef testVolumeSphereHexahedra < testUnfittedIntegration...
                                  & testUnfittedVolumeIntegration
    
   properties (Access = protected)
        testName = 'test_sphere_hexahedra';  
        analiticalArea = (4/3)*pi;
        meshType = 'INTERIOR';
   end
 
end

