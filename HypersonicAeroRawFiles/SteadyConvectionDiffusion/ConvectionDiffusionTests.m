classdef ConvectionDiffusionTests < handle & matlab.unittest.TestCase

    properties (TestParameter)
        tests  = {'lowPecletLinearElem','highPecletQuadElem'};
        p2test = {'highPecletQuadElem'};
    end

    methods (Test, TestTags = {'ConvectionDiffusion'})
        function testWithoutSourceTermGalerkin(testCase, tests)
            run(tests)
            sol = computeSteadyConvectionDiffusion1D(1,a,nu,100,p,1);
            load([tests,'Steady1DGalerkinNoSource.mat'],'x');
            err = norm(sol-x)/norm(x);
            tol = 1e-6;
            testCase.verifyLessThanOrEqual(err, tol)
            close all;
        end

        function testWithSourceTermGalerkin(testCase,tests)
            run(tests)
            sol = computeSteadyConvectionDiffusion1D(3,a,nu,100,p,1);
            load([tests,'Steady1DGalerkinSource.mat'],'x');
            err = norm(sol-x)/norm(x);
            tol = 1e-6;
            testCase.verifyLessThanOrEqual(err, tol)
            close all;
        end

        function testWithoutSourceTermSU(testCase,tests)
            run(tests)
            sol = computeSteadyConvectionDiffusion1D(1,a,nu,100,p,2);
            load([tests,'Steady1DSUNoSource.mat'],'x');
            err = norm(sol-x)/norm(x);
            tol = 1e-6;
            testCase.verifyLessThanOrEqual(err, tol)
            close all;
        end

        function testWithSourceTermSU(testCase,tests)
            run(tests)
            sol = computeSteadyConvectionDiffusion1D(3,a,nu,100,p,2);
            load([tests,'Steady1DSUSource.mat'],'x');
            err = norm(sol-x)/norm(x);
            tol = 1e-6;
            testCase.verifyLessThanOrEqual(err, tol)
            close all;
        end

        function testWithoutSourceTermSUPG(testCase,tests)
            run(tests)
            sol = computeSteadyConvectionDiffusion1D(1,a,nu,100,p,3);
            load([tests,'Steady1DSUPGNoSource.mat'],'x');
            err = norm(sol-x)/norm(x);
            tol = 1e-6;
            testCase.verifyLessThanOrEqual(err, tol)
            close all;
        end

        function testWithSourceTermSUPG(testCase,tests)
            run(tests)
            sol = computeSteadyConvectionDiffusion1D(3,a,nu,100,p,3);
            load([tests,'Steady1DSUPGSource.mat'],'x');
            err = norm(sol-x)/norm(x);
            tol = 1e-6;
            testCase.verifyLessThanOrEqual(err, tol)
            close all;
        end

        function testWithoutSourceTermGLS(testCase,p2test)
            run(p2test)
            sol = computeSteadyConvectionDiffusion1D(1,a,nu,100,p,4);
            load([p2test,'Steady1DGLSNoSource.mat'],'x');
            err = norm(sol-x)/norm(x);
            tol = 1e-6;
            testCase.verifyLessThanOrEqual(err, tol)
            close all;
        end

        function testWithSourceTermGLS(testCase,p2test)
            run(p2test)
            sol = computeSteadyConvectionDiffusion1D(3,a,nu,100,p,4);
            load([p2test,'Steady1DGLSSource.mat'],'x');
            err = norm(sol-x)/norm(x);
            tol = 1e-6;
            testCase.verifyLessThanOrEqual(err, tol)
            close all;
        end

        function testWithoutSourceTermSGS(testCase,p2test)
            run(p2test)
            sol = computeSteadyConvectionDiffusion1D(1,a,nu,100,p,5);
            load([p2test,'Steady1DSGSNoSource.mat'],'x');
            err = norm(sol-x)/norm(x);
            tol = 1e-6;
            testCase.verifyLessThanOrEqual(err, tol)
            close all;
        end

        function testWithSourceTermSGS(testCase,p2test)
            run(p2test)
            sol = computeSteadyConvectionDiffusion1D(3,a,nu,100,p,5);
            load([p2test,'Steady1DSGSource.mat'],'x');
            err = norm(sol-x)/norm(x);
            tol = 1e-6;
            testCase.verifyLessThanOrEqual(err, tol)
            close all;
        end

    end
end