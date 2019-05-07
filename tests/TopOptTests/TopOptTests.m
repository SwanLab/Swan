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
                'testGripping';
                'testMicro';
                'testCantilever';
                'testMicro2';
                'testProjectedSlerp';
                'testCantilever3';
                'testBridge2';
                'testCantilever2';
                'testBridge';
                
                };
        end
        
    end
    
end
