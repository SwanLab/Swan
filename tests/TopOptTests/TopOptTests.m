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
               'testCantilever3';
               'testProjectedSlerp';                
               'testBridge2';                   
               'testBridge';                  
               'testMicro';               
               'testCantilever2';                
               'testCantilever';                
               'testGripping';  
               'testMicro2';
               };

        end
    end
    
end
