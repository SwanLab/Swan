classdef DiffReactTests < matlab.unittest.TestCase

    properties (TestParameter)
%         file = {'testDiffReactHexagon', 'testDiffReactTorus', 'testDiffReactCube'}
        file = {'testDiffReactHexagon'}
        file3d = {'testDiffReactTorus', 'testDiffReactCube'}
        LHStype = {'StiffnessMass', 'StiffnessMassBoundaryMass'}
    end

    methods (Test, TestTags = {'DiffReact', '2D'})

        function testHexagon(testCase, file, LHStype)
            s   = testCase.createFEMparameters(file, LHStype);
            RHS = testCase.createRHS(s.mesh);
            fem = PhysicalProblem.create(s);
            fem.computeLHS(0.1857);
            fem.computeVariables(RHS);
%             fem.print(filename)
            err = testCase.computeError(file, LHStype, fem);
            tol = 1e-6;
            testCase.verifyLessThanOrEqual(err, tol)
        end

    end

    methods (Test, TestTags = {'DiffReact', '3D'})

        function test3D(testCase, file3d, LHStype)
            s   = testCase.createFEMparameters(file3d, LHStype);
            RHS = testCase.createRHS(s.mesh);
            fem = PhysicalProblem.create(s);
            fem.computeLHS(0.1857);
            fem.computeVariables(RHS);
%             fem.print(filename)
            err = testCase.computeError(file3d, LHStype, fem);
            tol = 1e-6;
            testCase.verifyLessThanOrEqual(err, tol)
        end

    end

    methods (Access = private)

        function error = computeError(testCase, file, LHStype, fem)
            file2load = append(file, '_', LHStype);
            cV = fem.variables.x;
            sV = load(file2load).x;
            error = norm(sV - cV)/norm(sV);
        end

        function name = loadFile(testCase, file, LHStype)
            name = append(file, '_', LHStype);
        end

        function s = createFEMparameters(testCase, file, LHStype)
            gidParams = testCase.createGiDparameters(file);
            s.dim       = gidParams.pdim;
            s.type      = gidParams.ptype;
            s.scale     = gidParams.scale;
            s.mesh      = gidParams.mesh;
            s.LHStype   = LHStype;
        end
        
        function gidParams = createGiDparameters(testCase, file)
            gidReader = FemInputReaderGiD();
            gidParams = gidReader.read(file);
        end
        
        function rhs = createRHS(testCase, mesh)
            M = testCase.computeM(mesh);
            u = testCase.createDisplacement(M);
            rhs = M*u;
        end
        
        function M = computeM(testCase, mesh)
            a.mesh    = mesh;
            a.fValues = zeros(mesh.nnodes, 1);
            a.order   = 'P1';
            f = LagrangianFunction(a);
            s.type         = 'MassMatrix';
            s.quadratureOrder = 2;
            s.mesh = mesh;
            s.test  = f;
            s.trial = f;
            LHS = LHSIntegrator.create(s);
            M = LHS.compute();
        end
        
        function u = createDisplacement(testCase, M)
            sizeM = size(M,1);
            sizeI = floor(sizeM/2);
            sizeJ = sizeM - sizeI;
            u = [ones(sizeI, 1); zeros(sizeJ, 1)];
        end

    end

end