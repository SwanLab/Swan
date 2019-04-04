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
                
            'testProjectedSlerp';
    
            'testMicro';
            'testCantilever2';            
            
            'testCantilever3';
                
            'testCantilever';
            'testGripping';
                
            'testBridge';
            'testMicro2';
            'testBridge2';
            

            
            };
        end
        
    end
    
end
