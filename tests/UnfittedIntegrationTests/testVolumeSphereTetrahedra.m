classdef testVolumeSphereTetrahedra < testUnfittedIntegration...
                                  & testUnfittedVolumeIntegration
    
   properties (Access = protected)
        testName = 'test_sphere_tetrahedra';  
        analiticalArea = (4/3)*pi;
        meshType = 'INTERIOR';
   end
end

