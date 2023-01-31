classdef DehomogenizingTests < handle & matlab.unittest.TestCase

    properties (TestParameter)
    end

    methods (Test, TestTags = {'Singularities'})

        function testDehomogenizingSingularities(testCase)
            s.testName = 'test_dehomogenizingSingularities';            
            test = DehomogenizingSingularitiesTest(s);
            passed = test.hasPassed();
            verifyTrue(testCase, passed)
        end

    end

    properties (Access = protected)
        FieldOfStudy = 'Dehomogenizing'
        tests
    end

end