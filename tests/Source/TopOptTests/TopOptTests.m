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
                'testDualNestedInPrimalWithProjectedGradient';  
                'testDualNestedInPrimalWithSlerp';                                
                'testStressM1M2';
                'testM1M2';
                'testCantilever2';
                'testCantilever';
                'testCantilever3';
                'testBridge';
                'testGripping';
                'testMicro';
                'testBridge2'; 
                'testMicro2';                                
                };
        end

    end

end
