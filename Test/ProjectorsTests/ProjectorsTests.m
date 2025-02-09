classdef ProjectorsTests < handle & matlab.unittest.TestCase

    properties (TestParameter)
        FEMtests = {
            'test_P0toP1','test_P1toP0','test_P0toP1disc','test_P1toP1disc',...
            'test_P0toP0','test_P1toP1'
            }
    end
        
    methods (Test, TestTags = {'Projectors', 'Various', 'FEMTests'})

        function testProjectorsInFEM(testCase, FEMtests)
            s.computerType     = 'PROJECT';
            s.problemType      = 'FEM';
            s.testName         = FEMtests;
            s.variablesToStore = {'xP'};
            test = PrecomputedVariableTest(s);
            err = test.computeError();
            tol = 1e-6;
            testCase.verifyLessThanOrEqual(err, tol)
        end

    end

end