classdef BoundaryCondTests < handle & matlab.unittest.TestCase

    properties (TestParameter)
        triangle = {'test2d_triangle'}
%         micro = {'test2d_micro', 'test2d_micro_thin',...
%             'test_micro_holeinclusion'}
        micro ={'test_micro_holeinclusion'}
    end

    methods (Test, TestTags = {'Monolithic', 'Macro'})

        function testTriangleNullDispMon(testCase, triangle)
            s.computerType     = 'FEM';
            s.testName         = triangle;
            s.variablesToStore = {'d_u'};
            s.computerType = 'FEM';
            s.solverType   = 'MONOLITHIC';
            s.solverMode   = 'DISP';
            test = PrecomputedVariableTest(s);
            err = test.computeError();
            tol = 1e-12;
            testCase.verifyLessThanOrEqual(err, tol);
        end

        function testTriangleNullDispRed(testCase, triangle)
            s.computerType     = 'FEM';
            s.testName         = triangle;
            s.variablesToStore = {'d_u'};
            s.computerType = 'FEM';
            s.solverType   = 'MONOLITHIC';
            s.solverMode   = 'DISP';
            test = PrecomputedVariableTest(s);
            err = test.computeError();
            tol = 1e-12;
            testCase.verifyLessThanOrEqual(err, tol);
        end

        function testTriangleDispMon(testCase, triangle)
            s.computerType     = 'FEM';
            s.testName         = [triangle '_non_null'];
            s.variablesToStore = {'d_u'};
            s.computerType = 'FEM';
            s.solverType   = 'MONOLITHIC';
            s.solverMode   = 'DISP';
            s.testResultsName  = [triangle '_non_null'];
            test = PrecomputedVariableTest(s);
            err = test.computeError();
            tol = 1e-12;
            testCase.verifyLessThanOrEqual(err, tol);
        end

        function testTriangleDispRed(testCase, triangle)
            s.computerType     = 'FEM';
            s.testName         = [triangle '_non_null'];
            s.variablesToStore = {'d_u'};
            s.computerType = 'FEM';
            s.solverType   = 'REDUCED';
            s.solverMode   = 'DISP';
            s.testResultsName  = [triangle '_non_null'];
            test = PrecomputedVariableTest(s);
            err = test.computeError();
            tol = 1e-12;
            testCase.verifyLessThanOrEqual(err, tol);
        end

    end

    methods (Test, TestTags = {'Monolithic', 'Micro'})

        function testMicroDispMonolitic(testCase, micro)
            s.testName = micro;
            s.variablesToStore = {'Chomog'};
            s.computerType = 'MICRO';
            s.solverType   = 'MONOLITHIC';
            s.solverMode   = 'DISP';
            test = PrecomputedVariableTest(s);
            err = test.computeError();
            tol = 1e-4;
            testCase.verifyLessThanOrEqual(err, tol)
        end

        function testMicroFlucReduced(testCase, micro)
            s.testName = micro;
            s.variablesToStore = {'Chomog'};
            s.computerType = 'MICRO';
            s.solverType   = 'REDUCED';
            s.solverMode   = 'FLUC';
            test = PrecomputedVariableTest(s);
            err = test.computeError();
            tol = 1e-12;
            testCase.verifyLessThanOrEqual(err, tol)
        end
        function testMicroFlucMonolitic(testCase, micro)
            s.testName = micro;
            s.variablesToStore = {'Chomog'};
            s.computerType = 'MICRO';
            s.solverType   = 'MONOLITHIC';
            s.solverMode   = 'FLUC';
            test = PrecomputedVariableTest(s);
            err = test.computeError();
            tol = 1e-12;
            testCase.verifyLessThanOrEqual(err, tol)
        end

    end
end