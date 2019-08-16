classdef UnfittedIntegrationTests < testRunner
    
    properties (Access = protected)
        FieldOfStudy = 'Unfitted integration'
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
                
            'testPerimeterRectangleTriangle'
            'testPerimeterRectangleQuadrilateral'
            
            'testAreaCircleQuadrilateral'
            'testAreaCircleTriangle'
            
            'testPerimeterCircleTriangle'
            'testPerimeterCircleQuadrilateral'
            
            'testSurfaceSphereTetrahedra';
            'testSurfaceSphereHexahedra';
            
            'testVolumeSphereTetrahedra';
            'testVolumeSphereHexahedra';
            
            'testSurfaceCylinderTetrahedra';
            'testSurfaceCylinderHexahedra';
            
            'testVolumeCylinderTetrahedra';
            'testVolumeCylinderHexahedra';
            };
        end
        
    end
    
end

