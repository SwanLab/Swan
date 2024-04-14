classdef ConvectionDiffusionTests < handle & matlab.unittest.TestCase

    properties (TestParameter)
        tests = {'lowPecletLinearElem'};
    end

    methods (Test, TestTags = {'ConvectionDiffusion','Steady1DGalerkinNoSource', 'Classic'})
        function testWithoutSourceTermGalerkin(testCase, tests)
            run(tests)
            sol = computeSteadyConvectionDiffusion1D(1,a,nu,100,p,1);
            load([tests,'Steady1DGalerkinNoSource.mat'],'x');
            err = norm(sol-x)/norm(x);
            tol = 1e-6;
            testCase.verifyLessThanOrEqual(err, tol)
            close all;
        end

    end
end