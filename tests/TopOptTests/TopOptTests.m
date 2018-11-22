
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
               'testCantilever';   
               'testMicro';
               'testMicro2';
               'testBridge2';   
               'testBridge';
               'testCantilever2';
               'testCantilever3';
               'testProjectedSlerp';
               'testGripping';              
               };

        end
    end
    
end

