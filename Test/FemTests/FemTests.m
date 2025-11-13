classdef FemTests < handle & matlab.unittest.TestCase

    properties (TestParameter)
        triangle = {'test2d_triangle'}
        quad = {'test2d_quad'}
        tests2d = {'test2d_triangle', 'test2d_quad'}
        tests3d = {'test3d_tetrahedra', 'test3d_hexahedra'}
        hexahedra = {'test3d_hexahedra'}
        duTests = {'test2d_triangle', 'test2d_quad', 'test3d_tetrahedra', 'test3d_hexahedra'}
        stokesTests = {'test2d_stokes_triangle_steady', 'test2d_stokes_triangle_transient'}
        microTests = {'test2d_micro', 'test3d_micro_cube'}
        thermalTests = {'test_thermal'}
        hyperelasticTests = {'test_hyperelastic'}
        solvers = {'test_pyAMG'};
    end

    methods (Test, TestTags = {'FEM','Solvers'})

        function testIterativeSolvers(testCase, solvers)
            load('Por3D_LHS.mat','LHS');
            load('Por3D_RHS.mat','RHS');
            load([solvers,'.mat'],'xReal');
            run(solvers);
            solver = Solver.create(s);
            xNew = solver.solve(LHS,RHS);
            err = norm(xReal - xNew)/norm(xReal);
            tol = 1e-6;
            testCase.verifyLessThanOrEqual(err, tol);
        end

    end

    methods (Test, TestTags = {'Triangle'})

        function testTriangle(testCase, triangle)
            s.computerType     = 'FEM';
            s.testName         = triangle;
            s.variablesToStore = {'d_u'};
            test = PrecomputedVariableTest(s);
            err = test.computeError();
            tol = 1e-6;
            testCase.verifyLessThanOrEqual(err, tol)
        end

    end

    methods (Test, TestTags = {'FEM', 'Quadratic'})

        function testTriangleQuadratic(testCase, triangle)
            s.computerType     = 'FEM';
            s.testName         = triangle;
            s.testResultsName  = [triangle '_quadratic'];
            s.variablesToStore = {'d_u'};
            s.interpolationType = 'QUADRATIC';
            test = PrecomputedVariableTest(s);
            err = test.computeError();
            tol = 1e-6;
            testCase.verifyLessThanOrEqual(err, tol)
        end

    end
% 
%     methods (Test, TestTags = {'Quad'})
% 
%         function testQuad(testCase, quad)
%             s.computerType     = 'FEM';
%             s.testName         = quad;
%             s.variablesToStore = {'d_u'};
%             test = PrecomputedVariableTest(s);
%             err = test.computeError();
%             tol = 1e-6;
%             testCase.verifyLessThanOrEqual(err, tol)
%         end
% 
%     end
% 
% 
%     methods (Test, TestTags = {'2D'})
% 
%         function test2d(testCase, tests2d)
%             s.computerType    = 'FEM';
%             s.testName         = tests2d;
%             s.variablesToStore = {'d_u'};
%             test = PrecomputedVariableTest(s);
%             err = test.computeError();
%             tol = 1e-6;
%             testCase.verifyLessThanOrEqual(err, tol)
%         end
% 
%     end
% 
% 
%     methods (Test, TestTags = {'3D'})
% 
%         function test3d(testCase, tests3d)
%             s.computerType    = 'FEM';
%             s.testName         = tests3d;
%             s.variablesToStore = {'d_u'};
%             test = PrecomputedVariableTest(s);
%             err = test.computeError();
%             tol = 1e-6;
%             testCase.verifyLessThanOrEqual(err, tol)
%         end
% 
%     end
% 
%     methods (Test, TestTags = {'Hexahedra'})
% 
%         function testHexahedra(testCase, hexahedra)
%             s.computerType    = 'FEM';
%             s.testName         = hexahedra;
%             s.variablesToStore = {'d_u'};
%             test = PrecomputedVariableTest(s);
%             err = test.computeError();
%             tol = 1e-6;
%             testCase.verifyLessThanOrEqual(err, tol)
%         end
% 
%     end
% 
%     methods (Test, TestTags = {'FEM', 'Passed', 'Classic', 'Displacement', 'ToPass'})
% 
%         function testDisplacement(testCase, duTests)
%             s.computerType    = 'FEM'; %FEM
%             s.testName         = duTests;
%             s.variablesToStore = {'d_u'};
%             test = PrecomputedVariableTest(s);
%             err = test.computeError();
%             tol = 1e-6;
%             testCase.verifyLessThanOrEqual(err, tol)
%         end
% 
%     end
% 
%     methods(Test, TestTags = {'FEM', 'Passed', 'Classic', 'Stokes'})
% 
%         function testStokes(testCase, stokesTests)
%             s.computerType     = 'STOKES';
%             s.testName         = stokesTests;
%             s.variablesToStore = {'u','p'};
%             test = PrecomputedVariableTest(s);
%             err = test.computeError();
%             tol = 1e-6;
%             testCase.verifyLessThanOrEqual(err, tol)
%         end
% 
%     end
% 
%     methods(Test, TestTags = {'FEM', 'Classic', 'Micro'})
% 
%         function testMicro(testCase, microTests)
%             s.testName = microTests;
%             s.variablesToStore = {'Chomog'};
%             s.computerType = 'MICRO';
%             test = PrecomputedVariableTest(s);
%             err = test.computeError();
%             tol = 1e-6;
%             testCase.verifyLessThanOrEqual(err, tol)
%         end
% 
%     end
% 
%     methods(Test, TestTags = {'Thermal'})
% 
%         function testThermal(testCase, thermalTests)
%             s.testName = thermalTests;
%             s.variablesToStore = {'d_u'};
%             s.computerType = 'THERMAL';
%             test = PrecomputedVariableTest(s);
%             err = test.computeError();
%             tol = 1e-6;
%             testCase.verifyLessThanOrEqual(err, tol)
%         end
% 
%     end
% 
%     methods(Test, TestTags = {'Hyperelastic'})
% 
%         function testHyperelastic(testCase, hyperelasticTests)
%             s.testName = hyperelasticTests;
%             s.variablesToStore = {'d_u'};
%             s.computerType = 'FEM';
%             test = PrecomputedVariableTest(s);
%             err = test.computeError();
%             tol = 1e-6;
%             testCase.verifyLessThanOrEqual(err, tol)
%         end
% 
%     end

    methods(Test, TestTags = {'FEM', 'Classic', 'PrincipalDirection'})

        function testPrincipalDirection2D(testCase)
            s.stressDim = 3;
            s.pdim      = 2;
            s.nGaus     = 3;
            s.mesh = TriangleMesh(1,1,5,5);
            test = PrincipalDirectionTest(s);
            err = test.computeError();
            tol = 1e-12;
            testCase.verifyLessThanOrEqual(err, tol)
        end

        function testPrincipalDirection3D(testCase)
            s.stressDim = 6;
            s.pdim      = 3;
            s.nGaus     = 3;
            s.mesh = TetraMesh(1,1,1,5,5,5);            
            test = PrincipalDirectionTest(s);
            err = test.computeError();
            tol = 1e-12;
            testCase.verifyLessThanOrEqual(err, tol)
        end

    end

end