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
               'testCantilever3';                   
               'testCantilever2';
               'testBridge';                             
               'testCantilever';                
               'testGripping';  
                                
               'testProjectedSlerp';                                
               'testMicro2'; 
               'testBridge2';                
               };

        end
    end
    
end
