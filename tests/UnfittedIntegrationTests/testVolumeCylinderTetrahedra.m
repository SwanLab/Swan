classdef testVolumeCylinderTetrahedra < testUnfittedIntegration...
                                  & testUnfittedVolumeIntegration
    
   properties (Access = protected)
        testName = 'test_cylinder_tetrahedra';  
        analiticalArea = pi*2;
        meshType = 'INTERIOR';
   end
 
end

