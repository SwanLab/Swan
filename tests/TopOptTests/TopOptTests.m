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
            'testCantilever3';
            'testBridge';
            'testGripping';
            'testCantilever';
            'testMicro2';
            'testProjectedSlerp';
            'testCantilever2';
            'testMicro';

            };
        end

    end

end
