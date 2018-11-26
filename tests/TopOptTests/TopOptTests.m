
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
               'testBridge2';   
               'testBridge';
               'testCantilever2';
               'testProjectedSlerp';
               'testGripping';              
               'testMicro2';
               'testCantilever';   
               'testMicro';
               };

        end
    end
    
end

