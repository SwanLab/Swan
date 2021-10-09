classdef NewUnfittedIntegrationTests < handle & matlab.unittest.TestCase

    properties (TestParameter)
        perimeterCircleTests = {'test_circle_triangle', 'test_circle_quadrilateral'}
        stokesTests = {'test2d_stokes_triangle'}
        microTests = {'test2d_micro'}
    end

    methods (Test, TestTags = {'UnfittedIntegration', 'Perimeter', 'Circle'}) % my eyes are circles

        function testPerimeterCircle(testCase, perimeterCircleTests)
            s.solver           = 'FEM_SOLVER';
            s.testName         = perimeterCircleTests;
            s.variablesToStore = {'d_u'};
            analyticalValue = 2*pi;
            meshType = 'BOUNDARY';
            meshIncludeBoxContour = false;
            test = PrecomputedVariableTest(s);
            err = test.computeError();
            tol = 1e-6;
            testCase.verifyLessThanOrEqual(err, tol)
        end

    end

end

