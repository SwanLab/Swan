classdef test_all_micro < testRunner
    
    
    properties (Access = protected)
        FieldOfStudy = 'AmplificatorTest'
        tests
    end
    
    methods (Access = public)
        function obj = test_all_micro()
            obj@testRunner();
        end
    end
    
    methods (Access = protected)
        function loadTests(obj)
            obj.tests = {...                
                  'test2d_micro';
                  'test_micro';
                  'test_micro2';
                };

        end
    end
    
end

