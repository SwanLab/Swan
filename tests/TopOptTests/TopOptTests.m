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
                
            'testBridge';
            'testMicro';
            'testCantilever3';
%             'testCantilever2'; ipopt peta
            
            'testCantilever';
            'testGripping';
            
            'testProjectedSlerp';
            'testBridge2';
            'testMicro2';

            };
        end
        
    end
    
end
