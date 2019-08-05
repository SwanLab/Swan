classdef ImageProcessingTests < testRunner
        
    properties (Access = protected)
        FieldOfStudy = 'ImageProcessing'
        tests
    end
    
    methods (Access = public)
        function obj = ImageProcessingTests()
            obj@testRunner();
        end
    end
    
    methods (Access = protected)
        function loadTests(obj)
            obj.tests = {...  
                'testDenoisingEinstein';
                'testAcceleratedDenoisingEinstein';
                };

        end
    end
    
end