classdef ProjectorsTests < handle & matlab.unittest.TestCase

    properties (TestParameter)
        fastDisp = {
            'test_P0toP1','test_P1toP0'
            }
    end
        
    methods (Test, TestTags = {'Projectors', 'Various', 'FastDisp'})

        function testFastDisplacement(testCase, fastDisp)
            s.computerType     = 'PROJECT';
            s.problemType      = 'FEM';
            s.testName         = fastDisp;
            s.variablesToStore = {'xP'};
            test = PrecomputedVariableTest(s);
            err = test.computeError();
            tol = 1e-6;
            testCase.verifyLessThanOrEqual(err, tol)
        end

    end

end