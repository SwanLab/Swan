classdef ZZZTopOptTests < handle & matlab.unittest.TestCase

    properties (TestParameter)
        compTests = {'test_bridge', 'test_bridge2', ...
            'test_cantilever', 'test_cantilever2', 'test_cantilever3', ...
            'test_micro', 'test_micro2', ...
            'testDualNestedInPrimal_WithProjectedGradient', ...
            'testDualNestedInPrimal_WithSlerp', ...
            'test_interiorPerimeter'}
    end

    methods (Test, TestTags = {'TopOpt', 'Passed', 'Classic', 'Displacement', 'Slow'})

        function testDisplacement(testCase, didnt)
            cd ../../../ % NOOOO, al segon test torna a saltar enrere
            s.computerType    = 'TOPOPT';
            s.testName         = didnt;
            s.variablesToStore = {'x'};
            test = PrecomputedVariableTest(s);
            err = test.computeError();
            tol = 1e-6;
            testCase.verifyLessThanOrEqual(err, tol)
            cd ./tests/Source/TopOptTests/
        end

    end

    methods (Test, TestTags = {'TopOpt', 'Passed', 'Classic', 'Displacement', 'Slow'})

        function testAnalyticVsRegularizedPerimeter(testCase, didnt)
            cd ../../../ % NOOOO, al segon test torna a saltar enrere
            s.computerType    = 'TOPOPT';
            s.testName         = didnt;
            s.variablesToStore = {'x'};
            test = PrecomputedVariableTest(s);
            err = test.computeError();
            tol = 1e-6;
            testCase.verifyLessThanOrEqual(err, tol)
            cd ./tests/Source/TopOptTests/
        end

    end

end

