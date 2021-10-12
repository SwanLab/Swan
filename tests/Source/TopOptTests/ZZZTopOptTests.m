classdef ZZZTopOptTests < handle & matlab.unittest.TestCase

    properties (TestParameter)
        compTestsNO = {'test_M1M2', 'test_StressM1M2'}
        compTests = {'test_bridge', 'test_bridge2', ...
            'test_cantilever', 'test_cantilever2', 'test_cantilever3', ...
            'test_micro', 'test_micro2', ...
            'testDualNestedInPrimal_WithProjectedGradient', ...
            'testDualNestedInPrimal_WithSlerp', ...
            'test_interiorPerimeter'}
        dimensions = {'2D', '3D'}
        vigdergauzTests = {'test_VigergauzMicroStructure', 'test_VigergauzMicroStructureFromStrain'}
        vigdergauzVolumes = {0.6, 0.75}
    end

    methods (Test, TestTags = {'TopOpt', 'Passed', 'Classic', 'Displacement', 'Slow'})

        function testDisplacement(testCase, compTests)
            cd ../../../
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
            cd ../../../
            test = testAnalyticVsRegularizedPerimeter();
            err = test.computeError();
            tol = 5e-2;
            testCase.verifyLessThanOrEqual(err, tol)
            cd ./tests/Source/TopOptTests/
        end

    end

    methods (Test, TestTags = {'TopOpt', 'Passed', 'Classic', 'Displacement', 'Slow'})

        function testSuperEllipseExponent(testCase)
            cd ../../../
            test = testSuperEllipseExponent();
            err = test.computeError();
            tol = 1e-12;
            testCase.verifyLessThanOrEqual(err, tol)
            cd ./tests/Source/TopOptTests/
        end

    end

    methods (Test, TestTags = {'TopOpt', 'Passed', 'Classic', 'Displacement', 'Slow'})

        function testSimp(testCase, dimensions)
            cd ../../../
            test = SimplAllTestExplicitVsImplicit(dimensions);
            err = test.computeError();
            tol = 1e-12;
            testCase.verifyLessThanOrEqual(err, tol)
            cd ./tests/Source/TopOptTests/
        end

    end

    methods (Test, TestTags = {'TopOpt', 'Passed', 'Nou', 'Displacement', 'Slow'})

        function testVigdergauz(testCase, vigdergauzTests, vigdergauzVolumes)
            cd ../../../
            s.testName = vigdergauzTests;
            s.volume = vigdergauzVolumes;
            test = ZZZtestVigdergauzMicroStructure(s);
            err = test.computeError();
            tol = 2*1e-1;
            testCase.verifyLessThanOrEqual(err, tol)
            cd ./tests/Source/TopOptTests/
        end

    end

end

