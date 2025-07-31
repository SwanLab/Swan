classdef HyperelasticityTests < handle & matlab.unittest.TestCase

    properties (TestParameter)
        cases = {'HoleDirich','Metamaterial'}
    end

    methods (Test, TestTags = {'DisplReact'})
        function test2D(testCase, cases)
            filename = ['testHyperElas',cases,'2D'];
            s.fileName = filename;
            s.nsteps = testCase.computeNumberOfSteps(cases);
            s.printing = false;
            s.bcCase = cases;
            s.meshGen = cases;
            h = HyperelasticProblem(s);
            uNew = h.uFun.fValues(:);
            rNew = h.rFun.fValues(:);
            load(filename,'uRef','rRef');
            errU = norm(uNew(:)-uRef(:))/norm(uRef(:));
            errR = norm(rNew(:)-rRef(:))/norm(rRef(:));
            err  = max(errU,errR);
            tol  = 1e-6;
            testCase.verifyLessThanOrEqual(err, tol)
        end
    end

    methods (Static, Access = private)

        function nsteps = computeNumberOfSteps(cases)
            switch cases
                case 'HoleDirich'
                    nsteps = 20;
                case 'Metamaterial'
                    nsteps = 75;
            end
        end

    end

end