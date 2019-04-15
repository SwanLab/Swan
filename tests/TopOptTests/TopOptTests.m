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
            'testMicro2';
            'testBridge2';
            
            'testProjectedSlerp';
            
            'testCantilever3';
            
            'testMicro';
            'testCantilever2';            
                
            'testCantilever';
            'testGripping';
            
            };
        end
        
    end
    
end
