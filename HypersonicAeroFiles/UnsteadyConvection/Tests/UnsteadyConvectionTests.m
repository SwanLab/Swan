classdef UnsteadyConvectionTests < handle & matlab.unittest.TestCase

    properties (TestParameter)
        tests  = {1,2,3,4,5,6,7};
    end

    methods (Test, TestTags = {'UnsteadyConvection'})
        function testDiffMeth(testCase, tests)
            s.method      = tests;
            s.revolutions = 1;
            s.timeSteps   = 120;
            prob          = UnsteadyConvectionProblem(s);
            u             = prob.compute();
            load(['UnsteadyConvection',char(string(tests)),'.mat'],'x');
            err = norm(u(:)-x(:))/norm(x(:));
            tol = 1e-6;
            testCase.verifyLessThanOrEqual(err, tol)
            close all;
        end
    end
end