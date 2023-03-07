classdef FunctionTests < handle & matlab.unittest.TestCase

    properties (TestParameter)
            passedTests = {...
                'test1DLHS';
                }
    end

%     methods (Test, TestTags = {'FunctionTests', 'ShowingError'})
% 
%         function testsError(testCase, errorTests)
%             testCase.fixFolder();
%             test = eval(errorTests);
%             err = test.computeError();
%             tol = test.tol;
%             testCase.verifyLessThanOrEqual(err, tol)
%         end
% 
%     end

    methods (Test, TestTags = {'FunctionTests'})

        function testsPassed(testCase, passedTests)
            test = eval(passedTests);
            passed = test.hasPassed();
            verifyTrue(testCase, passed)
        end

    end

end