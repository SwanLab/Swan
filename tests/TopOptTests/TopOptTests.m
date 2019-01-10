
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
               'testCantilever2';                
               'testProjectedSlerp';
               'testCantilever';                
               'testGripping';  
               'testMicro';               
               'testMicro2';
               'testCantilever3';
               'testBridge2';   
               'testBridge';   
               };

        end
    end
    
end

