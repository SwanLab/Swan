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
                'testSphereTetrahedra';
                'testSphereHexahedra';
                };
            
        end
    end
    
end

