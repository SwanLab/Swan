classdef MultimaterialTests < handle & matlab.unittest.TestCase

    properties (TestParameter)

    end

    methods (Test, TestTags = {'MultiMat'})
        function test3Materials(testCase)
            problem = MultimaterialTesting();
            ls      = problem.solve();
            load('resultMultiMat.mat','lsReal');
            err = norm(ls-lsReal)/norm(lsReal);
            tol      = 1e-6;
            testCase.verifyLessThanOrEqual(err, tol)
        end

    end
  
end