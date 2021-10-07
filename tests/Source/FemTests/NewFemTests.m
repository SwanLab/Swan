classdef NewFemTests < handle % nota: extends TestRunner de normal
    % https://www.mathworks.com/help/matlab/ref/matlab.unittest.constraints-package.html
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

        function runTests(obj)
            nTests = size(obj.tests,1);
            for itest = 1: nTests
                i = obj.tests{itest};
                test = eval(i);
                test.checkTestPassed(i);
            end
        end
    end
    
end

