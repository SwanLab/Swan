classdef ConvectionDiffusionTests < handle & matlab.unittest.TestCase

    properties (TestParameter)
        tests  = {'lowPecletLinearElem','highPecletQuadElem'};
        p1test = {'lowPecletLinearElem'};
    end

    methods (Test, TestTags = {'ConvectionDiffusion'})
        function testWithoutSourceTermGalerkin(testCase, tests)
            run(tests)
            s.problem = 1;
            s.numel   = 100;
            s.p       = p;
            s.stab    = 1;
            prob      = SteadyConvectionDiffusionProblem(s);
            sol       = prob.compute(a,nu);
            load([tests,'Steady1DGalerkinNoSource.mat'],'x');
            err = norm(sol-x)/norm(x);
            tol = 1e-6;
            testCase.verifyLessThanOrEqual(err, tol)
            close all;
        end

        function testWithSourceTermGalerkin(testCase,tests)
            run(tests)
            s.problem = 3;
            s.numel   = 100;
            s.p       = p;
            s.stab    = 1;
            prob      = SteadyConvectionDiffusionProblem(s);
            sol       = prob.compute(a,nu);
            load([tests,'Steady1DGalerkinSource.mat'],'x');
            err = norm(sol-x)/norm(x);
            tol = 1e-6;
            testCase.verifyLessThanOrEqual(err, tol)
            close all;
        end

        function testWithoutSourceTermSU(testCase,tests)
            run(tests)
            s.problem = 1;
            s.numel   = 100;
            s.p       = p;
            s.stab    = 2;
            prob      = SteadyConvectionDiffusionProblem(s);
            sol       = prob.compute(a,nu);
            load([tests,'Steady1DSUNoSource.mat'],'x');
            err = norm(sol-x)/norm(x);
            tol = 1e-6;
            testCase.verifyLessThanOrEqual(err, tol)
            close all;
        end

        function testWithSourceTermSU(testCase,tests)
            run(tests)
            s.problem = 3;
            s.numel   = 100;
            s.p       = p;
            s.stab    = 2;
            prob      = SteadyConvectionDiffusionProblem(s);
            sol       = prob.compute(a,nu);
            load([tests,'Steady1DSUSource.mat'],'x');
            err = norm(sol-x)/norm(x);
            tol = 1e-6;
            testCase.verifyLessThanOrEqual(err, tol)
            close all;
        end

        function testWithoutSourceTermSUPG(testCase,p1test)
            run(p1test)
            s.problem = 1;
            s.numel   = 100;
            s.p       = p;
            s.stab    = 3;
            prob      = SteadyConvectionDiffusionProblem(s);
            sol       = prob.compute(a,nu);
            load([p1test,'Steady1DSUPGNoSource.mat'],'x');
            err = norm(sol-x)/norm(x);
            tol = 1e-6;
            testCase.verifyLessThanOrEqual(err, tol)
            close all;
        end

        function testWithSourceTermSUPG(testCase,p1test)
            run(p1test)
            s.problem = 3;
            s.numel   = 100;
            s.p       = p;
            s.stab    = 3;
            prob      = SteadyConvectionDiffusionProblem(s);
            sol       = prob.compute(a,nu);
            load([p1test,'Steady1DSUPGSource.mat'],'x');
            err = norm(sol-x)/norm(x);
            tol = 1e-6;
            testCase.verifyLessThanOrEqual(err, tol)
            close all;
        end

    end
end