classdef (TestTags = {'FEM'}) ...
        NewFemTests < handle & matlab.unittest.TestCase
    % https://www.mathworks.com/help/matlab/ref/matlab.unittest.constraints-package.html
    properties (TestParameter)
        paramtests = {'NewMatlabTest2dTriangle', 'NewMatlabTest3dTetrahedra'};
        noms = {'test2d_triangle', 'test3d_tetrahedra'}
    end

    properties (Access = protected)
        FieldOfStudy = 'FEM'
        tests
    end
    
    methods (Test, TestTags = {'Passed', 'provesvelles'})
        function testPassed(testCase, paramtests)
            inst = eval(paramtests);
            err = inst.computeErrorForTest();
            tol = 1e-6;
            testCase.verifyLessThanOrEqual(err, tol)
        end
    end 

    methods (Test, TestTags = {'Passed', 'Legacy'})
        function testPassedNou(testCase, noms)
            s.testName = noms;
            s.variablesToStore = {'d_u'};
            inst = NewMatlabTest(s);
            err = inst.computeErrorForTest();
            tol = 1e-6;
            testCase.verifyLessThanOrEqual(err, tol)
        end
    end 

end

