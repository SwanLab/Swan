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
               'testBridge';
               'testMicro';
               
               'testBridge2';                
               'testCantilever2';
               'testProjectedSlerp';
               'testGripping';              
               'testMicro2';
               'testCantilever';   
               };
        end
    end  
end

