classdef ZZZTopOptTests < handle & matlab.unittest.TestCase

    properties (TestParameter)
        compTestsNO = {'test_M1M2', 'test_StressM1M2'}
        compTests = {'test_bridge', 'test_bridge2', ...
            'test_cantilever', 'test_cantilever2', 'test_cantilever3', ...
            'test_micro', 'test_micro2', ...
            'testDualNestedInPrimal_WithProjectedGradient', ...
            'testDualNestedInPrimal_WithSlerp', ...
            'test_interiorPerimeter'}
        explicitImplicitTests = {'2D', '3D'}
    end

    methods (Test, TestTags = {'TopOpt', 'Passed', 'Classic', 'Displacement', 'Slow'})

        function testDisplacement(testCase, compTests)
            cd ../../../ % NOOOO, al segon test torna a saltar enrere
            s.computerType    = 'TOPOPT';
            s.testName         = compTests;
            s.variablesToStore = {'x'};
            test = PrecomputedVariableTest(s);
            err = test.computeError();
            tol = 1e-6;
            testCase.verifyLessThanOrEqual(err, tol)
            cd ./tests/Source/TopOptTests/
        end

    end

    methods (Test, TestTags = {'TopOpt', 'Passed', 'Classic', 'Displacement', 'Slow'})

        function testAnalyticVsRegularizedPerimeter(testCase)
            cd ../../../ % NOOOO, al segon test torna a saltar enrere
            test = testAnalyticVsRegularizedPerimeter();
            err = test.computeError();
            tol = 5e-2;
            testCase.verifyLessThanOrEqual(err, tol)
            cd ./tests/Source/TopOptTests/
        end

    end

    methods (Test, TestTags = {'TopOpt', 'Passed', 'Classic', 'Displacement', 'Slow'})

        function testSuperEllipseExponent(testCase)
            cd ../../../ % NOOOO, al segon test torna a saltar enrere
            test = testSuperEllipseExponent();
            err = test.computeError();
            tol = 1e-12;
            testCase.verifyLessThanOrEqual(err, tol)
            cd ./tests/Source/TopOptTests/
        end

    end

    methods (Test, TestTags = {'TopOpt', 'Passed', 'Classic', 'Displacement', 'Slow'})

        function testSimp(testCase, explicitImplicitTests)
            cd ../../../ % NOOOO, al segon test torna a saltar enrere
            test = SimplAllTestExplicitVsImplicit(explicitImplicitTests);
            err = test.computeError();
            tol = 1e-12;
            testCase.verifyLessThanOrEqual(err, tol)
            cd ./tests/Source/TopOptTests/
        end

    end

end

