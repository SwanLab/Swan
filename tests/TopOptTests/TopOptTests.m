classdef TopOptTests < testRunner
    
    properties (Access = protected)
        FieldOfStudy = 'Topology Optimization'
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
                'testBridge';
                'testBridge2';
                'testDualNestedInPrimalWithSlerp';
                'testDualNestedInPrimalWithProjectedGradient';
                'testCantilever2';
                'testMicro2';
                'testCantilever3';
                'testGripping';
                'testMicro';
                'testCantilever';
                };
        end
        
    end
    
end
