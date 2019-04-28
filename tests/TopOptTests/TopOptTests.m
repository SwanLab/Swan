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
            'testMicro';                
            'testProjectedSlerp';                                                   
            'testCantilever';   
            'testCantilever2';                        
            'testGripping';                
            'testBridge2';                
            'testCantilever3';
            'testBridge';           
            'testMicro2';            
            };
        end

    end

end
