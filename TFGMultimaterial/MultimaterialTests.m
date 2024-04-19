classdef MultimaterialTests < handle & matlab.unittest.TestCase

    properties (TestParameter)

    end

    methods (Test, TestTags = {'MultiMat'})
        function test3Materials(testCase)
            psiC = mainTFGMultimaterial();
            load('resultMultiMat.mat','psiReal');
            err = norm(psiC-psiReal)/norm(psiReal);
            tol      = 1e-6;
            testCase.verifyLessThanOrEqual(err, tol)
        end

    end

end