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
                'testSurfaceCylinderTetrahedra';
                'testSurfaceCylinderHexahedra';
                
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
                
                'testVolumeCylinderTetrahedra';
                'testVolumeCylinderHexahedra';
                };
        end
        
    end
    
end

