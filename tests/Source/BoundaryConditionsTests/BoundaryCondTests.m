classdef BoundaryCondTests < handle & matlab.unittest.TestCase

    properties (TestParameter)
        triangle = {'test2d_triangle'}
%         micro = {'test2d_micro', 'test2d_micro_thin',...
%             'test_micro_holeinclusion'}
        micro ={'test_micro_holeinclusion'}
    end

    methods (Test, TestTags = {'Monolitic', 'Macro'})

        function testTriangleNullDispMon(testCase, triangle)
            s.computerType     = 'FEM';
            s.testName         = triangle;
            s.variablesToStore = {'d_u'};
            s.computerType = 'FEM';
            s.solType   = 'MONOLITIC';
            s.solMode   = 'DISP';
            test = PrecomputedVariableTest(s);
            err = test.computeError();
            tol = 1e-6;
            testCase.verifyLessThanOrEqual(err, tol);
        end

        function testTriangleNullDispRed(testCase, triangle)
            s.computerType     = 'FEM';
            s.testName         = triangle;
            s.variablesToStore = {'d_u'};
            s.computerType = 'FEM';
            s.solType   = 'REDUCED';
            s.solMode   = 'DISP';
            test = PrecomputedVariableTest(s);
            err = test.computeError();
            tol = 1e-6;
            testCase.verifyLessThanOrEqual(err, tol);
        end

        function testTriangleDispMon(testCase, triangle)
            s.computerType     = 'FEM';
            s.testName         = [triangle '_non_null'];
            s.variablesToStore = {'d_u'};
            s.computerType = 'FEM';
            s.solType   = 'MONOLITIC';
            s.solMode   = 'DISP';
            s.testResultsName  = [triangle '_non_null'];
            test = PrecomputedVariableTest(s);
            err = test.computeError();
            tol = 1e-6;
            testCase.verifyLessThanOrEqual(err, tol);
        end

        function testTriangleDispRed(testCase, triangle)
            s.computerType     = 'FEM';
            s.testName         = [triangle '_non_null2'];
            s.variablesToStore = {'d_u'};
            s.computerType = 'FEM';
            s.solType   = 'REDUCED';
            s.solMode   = 'DISP';
            s.testResultsName  = [triangle '_non_null2'];
            test = PrecomputedVariableTest(s);
            err = test.computeError();
            tol = 1e-6;
            testCase.verifyLessThanOrEqual(err, tol);
        end

        % function testTriangleDispMon_Original(testCase, triangle)
        %     s.computerType     = 'FEM';
        %     s.testName         = [triangle '_non_null'];
        %     s.variablesToStore = {'d_u'};
        %     s.computerType = 'FEM';
        %     s.solType   = 'MONOLITIC';
        %     s.solMode   = 'DISP';
        %     s.testResultsName  = [triangle '_non_null'];
        %     test = PrecomputedVariableTest(s);
        %     err = test.computeError();
        %     tol = 1e-6;
        %     testCase.verifyLessThanOrEqual(err, tol);
        % end
        % 
        % function testTriangleDispRed_Original(testCase, triangle)
        %     s.computerType     = 'FEM';
        %     s.testName         = [triangle '_non_null'];
        %     s.variablesToStore = {'d_u'};
        %     s.computerType = 'FEM';
        %     s.solType   = 'REDUCED';
        %     s.solMode   = 'DISP';
        %     s.testResultsName  = [triangle '_non_null'];
        %     test = PrecomputedVariableTest(s);
        %     err = test.computeError();
        %     tol = 1e-6;
        %     testCase.verifyLessThanOrEqual(err, tol);
        % end

    end

    methods (Test, TestTags = {'Monolitic', 'Micro'})

        function testMicroFlucReduced(testCase, micro)
            s.testName = micro;
            s.variablesToStore = {'Chomog'};
            s.computerType = 'MICRO';
            s.solType   = 'REDUCED';
            s.solMode   = 'FLUC';
            test = PrecomputedVariableTest(s);
            err = test.computeError();
            tol = 1e-4;
            testCase.verifyLessThanOrEqual(err, tol)
        end

        function testMicroDispMonolitic(testCase, micro)
            s.testName = micro;
            s.variablesToStore = {'Chomog'};
            s.computerType = 'MICRO';
            s.solType   = 'MONOLITIC';
            s.solMode   = 'DISP';
            test = PrecomputedVariableTest(s);
            err = test.computeError();
            tol = 1e-4;
            testCase.verifyLessThanOrEqual(err, tol)
        end

        function testMicroFlucMonolitic(testCase, micro)
            s.testName = micro;
            s.variablesToStore = {'Chomog'};
            s.computerType = 'MICRO';
            s.solType   = 'MONOLITIC';
            s.solMode   = 'FLUC';
            test = PrecomputedVariableTest(s);
            err = test.computeError();
            tol = 1e-4;
            testCase.verifyLessThanOrEqual(err, tol)
        end

    end
end