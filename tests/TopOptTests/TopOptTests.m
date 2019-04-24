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
            'testCantilever3';

            'testBridge2';

            'testGripping';

            'testCantilever';
            'testMicro2';

            'testProjectedSlerp';


            'testCantilever2';

            'testBridge';

            'testMicro';

            };
        end

    end

end
