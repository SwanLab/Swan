classdef MicroTests < testRunner
    
    
    properties (Access = protected)
        FieldOfStudy = 'AmplificatorTest'
        tests
    end
    
    methods (Access = public)
        function obj = MicroTests()
            obj@testRunner();
        end
    end
    
    methods (Access = protected)
        function loadTests(obj)
            obj.tests = {...                
                  'test2dMicro';
                  'testMicro';
                  'testMicro2';
                };

        end
    end
    
end

