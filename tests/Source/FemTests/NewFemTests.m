classdef (TestTags = {'FEM'}) ...
        NewFemTests < handle & matlab.unittest.TestCase
    % https://www.mathworks.com/help/matlab/ref/matlab.unittest.constraints-package.html
    properties (TestParameter)
        paramtests = {'NewMatlabTest2dTriangle', 'NewMatlabTest3dTetrahedra'};
    end

    properties (Access = protected)
        FieldOfStudy = 'FEM'
        tests
    end
    
    methods (Test, TestTags = {'Passed', 'Legacy'})
        function testPassed(testCase, paramtests) % 3 tests
            inst = eval(paramtests);
            err = inst.computeErrorForTest();
            tol = 1e-6;
            testCase.verifyLessThanOrEqual(err, tol)
        end
    end 
    
end

