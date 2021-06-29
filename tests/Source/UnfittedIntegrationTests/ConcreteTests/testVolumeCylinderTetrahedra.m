classdef testVolumeCylinderTetrahedra < testUnfittedIntegration_ExternalIntegrator...
                                  & testUnfittedVolumeIntegration
    
   properties (Access = protected)
        testName = 'test_cylinder_tetrahedra';  
        analyticalValue = pi*2;
        meshType = 'INTERIOR';
   end
 
end

