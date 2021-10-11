classdef FemTests < handle & matlab.unittest.TestCase

    properties (TestParameter)
        duTests = {'test2d_triangle', 'test2d_quad', 'test3d_tetrahedra', 'test3d_hexahedra'}
        stokesTests = {'test2d_stokes_triangle'}
        microTests = {'test2d_micro'}
    end

    methods (Test, TestTags = {'FEM', 'Passed', 'Classic', 'Displacement'})

        function testDisplacement(testCase, duTests)
            s.solver           = 'FEM_SOLVER';
            s.testName         = duTests;
            s.variablesToStore = {'d_u'};
            test = PrecomputedVariableTest(s);
            err = test.computeError();
            tol = 1e-6;
            testCase.verifyLessThanOrEqual(err, tol)
        end

    end

    methods(Test, TestTags = {'FEM', 'Passed', 'Classic', 'Stokes'})

        function testStokes(testCase, stokesTests)
            s.solver           = 'FEM_SOLVER';
            s.testName         = stokesTests;
            s.variablesToStore = {'u','p'};
            inst = StokesTest(s);
            err = inst.computeError();
            tol = 1e-6;
            testCase.verifyLessThanOrEqual(err, tol)
        end

    end

    methods(Test, TestTags = {'FEM', 'Classic', 'Micro'})

        function testMicro(testCase, microTests)
            s.testName = microTests;
            s.variablesToStore = {'Chomog'};
            s.solver = 'FEM_SOLVER';
            inst = Micro2DTest(s);
            err = inst.computeError();
            tol = 1e-6;
            testCase.verifyLessThanOrEqual(err, tol)
        end

    end

    methods(Test, TestTags = {'FEM', 'Classic', 'PrincipalDirection'})

        function testPrincipalDirection2D(testCase)
            s.stressDim = 3;
            s.pdim      = '2D';
            s.nelem     = 6400;
            s.nGaus     = 3;
            inst = PrincipalDirectionTest(s);
            err = inst.computeError();
            tol = 1e-12;
            testCase.verifyLessThanOrEqual(err, tol)
        end

        function testPrincipalDirection3D(testCase)
            s.stressDim = 6;
            s.pdim      = '3D';
            s.nelem     = 6400;
            s.nGaus     = 3;
            inst = PrincipalDirectionTest(s);
            err = inst.computeError();
            tol = 1e-12;
            testCase.verifyLessThanOrEqual(err, tol)
        end

    end

end

