classdef (TestTags = {'FEM'}) ...
        NewFemTests < handle & matlab.unittest.TestCase
    % https://www.mathworks.com/help/matlab/ref/matlab.unittest.constraints-package.html
    properties (TestParameter)
        paramtests = {'NewMatlabTest2dTriangle', 'NewMatlabTest3dTetrahedra'};
        duTests = {'test2d_triangle', 'test2d_quad', 'test3d_tetrahedra', 'test3d_hexahedra'}
        stokesTests = {'test2d_stokes_triangle'}
    end

    properties (Access = protected)
        FieldOfStudy = 'FEM'
        tests
    end
    
    methods (Test, TestTags = {'Passed', 'Legacy'})
        function testPassedLegacy(testCase, paramtests)
            inst = eval(paramtests);
            err = inst.computeErrorForTest();
            tol = 1e-6;
            testCase.verifyLessThanOrEqual(err, tol)
        end
    end 

    methods (Test, TestTags = {'Passed', 'Nou', 'du'})
        function testPassed(testCase, duTests)
            s.testName = duTests;
            s.variablesToStore = {'d_u'};
            s.solver = 'FEM_SOLVER';
            inst = NewMatlabTest(s);
            err = inst.computeErrorForTest();
            tol = 1e-6;
            testCase.verifyLessThanOrEqual(err, tol)
        end
    end 


    methods (Test, TestTags = {'Passed', 'Nou', 'stokes'})
        function testPassedStokes(testCase, stokesTests)
            s.testName = stokesTests;
%             s.variablesToStore = {'variable.u','variable.p'};
            s.variablesToStore = {'u','p'};
            s.solver = 'FEM_SOLVER';
            inst = NewMatlabTestStokes(s);
            err = inst.computeErrorForTest();
            tol = 1e-6;
            testCase.verifyLessThanOrEqual(err, tol)
        end
    end 
end

