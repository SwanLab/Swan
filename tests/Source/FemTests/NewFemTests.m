classdef NewFemTests < handle & matlab.unittest.TestCase % nota: extends TestRunner de normal
    % https://www.mathworks.com/help/matlab/ref/matlab.unittest.constraints-package.html
    properties (TestParameter)
        paramtests = {'NewMatlabTest2dTriangle'};
    end

    properties (Access = protected)
        FieldOfStudy = 'FEM'
        tests
    end
    
    methods (Access = public)
        function obj = NewFemTests()
            obj.performTest();
        end
    end
    
    methods (Access = private)
        function performTest(obj)
            obj.loadTests();
            obj.runTests();
        end
    end

    methods (Access = protected)

        function loadTests(obj)
            obj.tests = {'NewMatlabTest2dTriangle'};
        end

        function runSingleTest(obj, teststr)
            test = eval(teststr);
            test.checkTestPassed(teststr);
        end

        function runTests(obj)
            nTests = size(obj.tests,1);
            for itest = 1: nTests
                i = obj.tests{itest};
                test = eval(i);
                test.checkTestPassed(i);
            end
        end
    end
    
    methods (Test)
        function testNumel(testCase, paramtests) % 3 tests
            inst = eval(paramtests);
            err = inst.computeErrorForTest();
            tol = 1e-6;
            testCase.verifyLessThanOrEqual(err, tol)
        end
    end 
    
end

