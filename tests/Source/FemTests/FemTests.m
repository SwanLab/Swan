classdef FemTests < handle & matlab.unittest.TestCase

    properties (TestParameter)
        triangle = {'test2d_triangle'}
        quad = {'test2d_quad'}
        tests2d = {'test2d_triangle', 'test2d_quad'}
        tests3d = {'test3d_tetrahedra', 'test3d_hexahedra'}
        duTests = {'test2d_triangle', 'test2d_quad', 'test3d_tetrahedra', 'test3d_hexahedra'}
        stokesTests = {'test2d_stokes_triangle'}
        microTests = {'test2d_micro'}
    end

    methods (Test, TestTags = {'Triangle'})

        function testTriangle(testCase, triangle)
            s.computerType    = 'NEWFEM';
            s.testName         = triangle;
            s.variablesToStore = {'d_u'};
            test = PrecomputedVariableTest(s);
            err = test.computeError();
            tol = 1e-6;
            testCase.verifyLessThanOrEqual(err, tol)
        end

    end

    methods (Test, TestTags = {'Quad'})

        function testQuad(testCase, quad)
            s.computerType    = 'NEWFEM';
            s.testName         = quad;
            s.variablesToStore = {'d_u'};
            test = PrecomputedVariableTest(s);
            err = test.computeError();
            tol = 1e-6;
            testCase.verifyLessThanOrEqual(err, tol)
        end

    end


    methods (Test, TestTags = {'2D'})

        function test2d(testCase, tests2d)
            s.computerType    = 'NEWFEM';
            s.testName         = tests2d;
            s.variablesToStore = {'d_u'};
            test = PrecomputedVariableTest(s);
            err = test.computeError();
            tol = 1e-6;
            testCase.verifyLessThanOrEqual(err, tol)
        end

    end


    methods (Test, TestTags = {'3D'})

        function test3d(testCase, tests3d)
            s.computerType    = 'NEWFEM';
            s.testName         = tests3d;
            s.variablesToStore = {'d_u'};
            test = PrecomputedVariableTest(s);
            err = test.computeError();
            tol = 1e-6;
            testCase.verifyLessThanOrEqual(err, tol)
        end

    end
    methods (Test, TestTags = {'FEM', 'Passed', 'Classic', 'Displacement', 'ToPass'})

        function testDisplacement(testCase, duTests)
            s.computerType    = 'NEWFEM';
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
            s.computerType     = 'STOKES';
            s.testName         = stokesTests;
            s.variablesToStore = {'u','p'};
            test = PrecomputedVariableTest(s);
            err = test.computeError();
            tol = 1e-6;
            testCase.verifyLessThanOrEqual(err, tol)
        end

    end

    methods(Test, TestTags = {'FEM', 'Classic', 'Micro'})

        function testMicro(testCase, microTests)
            s.testName = microTests;
            s.variablesToStore = {'Chomog'};
            s.computerType = 'MICRO';
            test = PrecomputedVariableTest(s);
            err = test.computeError();
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
            test = PrincipalDirectionTest(s);
            err = test.computeError();
            tol = 1e-12;
            testCase.verifyLessThanOrEqual(err, tol)
        end

        function testPrincipalDirection3D(testCase)
            s.stressDim = 6;
            s.pdim      = '3D';
            s.nelem     = 6400;
            s.nGaus     = 3;
            test = PrincipalDirectionTest(s);
            err = test.computeError();
            tol = 1e-12;
            testCase.verifyLessThanOrEqual(err, tol)
        end

    end

end