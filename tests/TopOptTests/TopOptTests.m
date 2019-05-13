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
                
            'testMicro2';
            'testCantilever3';
            'testBridge2';
            'testBridge';
            'testGripping';
            'testMicro';
            'testCantilever';
            'testDualNestedInPrimalWithSlerp';
            'testDualNestedInPrimalWithProjectedGradient';
            'testCantilever2';
            };
        end
        
    end
    
end
