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
            'testMicro2';             
            'testCantilever2';                
            'testCantilever';                
            'testGripping';                
            'testBridge2';                
            'testCantilever3';
            'testBridge';
            'testProjectedSlerp';                           
            };
        end

    end

end
