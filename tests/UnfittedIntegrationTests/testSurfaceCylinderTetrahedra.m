classdef testSurfaceCylinderTetrahedra < testUnfittedIntegration...
                                  & testUnfittedSurfaceIntegration
    
   properties (Access = protected)
        testName = 'test_cylinder_tetrahedra';  
        analiticalArea = pi*2 + 2*pi*2;
        meshType = 'BOUNDARY';
   end
 
end

