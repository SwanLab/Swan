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
               'testMicro';                  
               'testProjectedSlerp';                                
               'testMicro2';               
               'testCantilever3';
               'testBridge2';                   
               'testBridge';                  
               'testCantilever2';                
               'testCantilever';                
               'testGripping';  
               };

        end
    end
    
end
