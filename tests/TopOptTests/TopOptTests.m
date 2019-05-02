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
            'testProjectedSlerp';                                                                                                                                                       
            'testBridge';                  
            'testGripping';
            'testCantilever2';                               
            'testMicro2';                      
            'testCantilever3';                
            'testBridge2';                                      
            'testMicro';                
            'testCantilever';         
            };
        end

    end

end
