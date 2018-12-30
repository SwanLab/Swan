classdef UnfittedIntegrationTests < testRunner
    
    properties (Access = protected)
        FieldOfStudy = 'Unfitted integration tests'
        tests
    end
    
    methods (Access = public)
        function  obj = UnfittedIntegrationTests()
            obj@testRunner();
        end
        
    end
    
    methods (Access = protected)
        function loadTests(obj)
            obj.tests = {...
                'testSurfaceCylinderTetrahedra';
                'testSurfaceCylinderHexahedra';
                
                'testVolumeCylinderTetrahedra';
                'testVolumeCylinderHexahedra';
                
                'testPerimeterCircleTriangle'  
                'testPerimeterCircleQuadrilateral'  
                
                'testAreaCircleTriangle'  
                'testAreaCircleQuadrilateral'  
                
                'testSurfaceSphereTetrahedra';
                'testSurfaceSphereHexahedra';
                
                'testVolumeSphereTetrahedra';
                'testVolumeSphereHexahedra';
                
                
                 
                };
            
%             'test_circle_triangle','test_circle_quadrilateral','test_sphere_tetrahedra','test_sphere_hexahedra','test_cylinder_tetrahedra','test_cylinder_hexahedra'};
        end
    end
    
end

