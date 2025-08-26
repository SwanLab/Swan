classdef HyperelasticityTests < handle & matlab.unittest.TestCase

    properties (TestParameter)
        cases = {'HoleDirich','Metamaterial'}
    end

    methods (Test, TestTags = {'DisplReact'})

        function test2D(testCase, cases)
            filename = ['testHyperElas',cases,'2D'];
            load(filename,'input','uRef','rRef');
            h = TestingHyperelasticity(input);
            output = h.compute();
            uNew = output.displacement.function{end}.fValues(:);
            rNew = output.reaction.function{end}.fValues(:);
            errU = norm(uNew(:)-uRef(:))/norm(uRef(:));
            errR = norm(rNew(:)-rRef(:))/norm(rRef(:));
            err  = max(errU,errR);
            tol  = 1e-6;
            testCase.verifyLessThanOrEqual(err, tol)
        end

    end

end