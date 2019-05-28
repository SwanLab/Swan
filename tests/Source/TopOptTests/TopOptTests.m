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
                'testBridge2';
                'testStressM1M2';
                'testM1M2';
                'testDualNestedInPrimalWithSlerp';
                'testDualNestedInPrimalWithProjectedGradient';
                'testMicro2';
                'testCantilever2';
                'testCantilever';
                'testCantilever3';
                'testBridge';
                'testGripping';
                'testMicro';
                };
        end

    end

end
