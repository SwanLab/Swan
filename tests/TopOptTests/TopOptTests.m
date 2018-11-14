classdef TopOptTests < testRunner
    
    
    properties (Access = protected)
        FieldOfStudy = 'TopOpt tests'
        tests
    end
    
    methods (Access = public)
        function obj = TopOptTests()
            obj@testRunner();
        end
    end
    
    methods (Access = protected)
        function loadTests(obj)
            obj.tests = {...    
               'testGripping';
               'testMicro';
               'testMicro2';
               'testBridge2';
               'testCantilever';
               'testCantilever2';
               'testCantilever3';
               'testProjectedSlerp';
               'testBridge';
                };

        end
    end
    
end

