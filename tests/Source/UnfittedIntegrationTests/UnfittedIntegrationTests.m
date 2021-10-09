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
            
            % agrupables analytical 2pi
           'MatlabtestPerimeterCircleTriangle'; % externalInt + perimeter
%            'testPerimeterCircleQuadrilateral'; % externalInt + perimeter
%            
%             % agrupables analytical pi
%            'testAreaCircleTriangle';
%            'testAreaCircleQuadrilateral';
%             
%            
%             % agrupables analytical 6
%            'testPerimeterRectangleTriangle';
%            'testPerimeterRectangleQuadrilateral';
%             
%             %%
%             % agrupables analytical 4pi
%            'testSurfaceSphereHexahedra';
%           'testSurfaceSphereTetrahedra';
%             % agrupables analytical 4/3 *pi
%            'testVolumeSphereTetrahedra';
%            'testVolumeSphereHexahedra';
%             
%              'testSurfaceCylinderTetrahedra';
%            'testSurfaceCylinderHexahedra';
%            'testVolumeCylinderTetrahedra';
%            'testVolumeCylinderHexahedra';
            
            
            };
        end
        
    end
    
end

