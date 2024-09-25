classdef ContinuumDamageComputer < handle

    properties (Access = private)
        mesh
        boundaryConditions
        material
    end

    methods (Access = public)

        function obj = ContinuumDamageComputer(cParams)
            obj.init(cParams)
        end

        function results = compute(obj)
            displacementFun = LagrangianFunction.create(obj.mesh, obj.mesh.ndim, 'P1');
            K = obj.computeK(displacementFun);
            F = obj.computeF(displacementFun);
            results = obj.computeU(K,F);
        end
    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh = cParams.mesh;
            obj.boundaryConditions = cParams.boundaryConditions;
            obj.material = cParams.material;
        end

        function K = computeK(obj,dispFun)

            ndimf = dispFun.ndimf;
            s.type     = 'ElasticStiffnessMatrix';
            s.mesh     = obj.mesh;
            s.test     = LagrangianFunction.create(obj.mesh,ndimf, 'P1');
            s.trial    = dispFun;
            s.material = obj.material;
            s.quadratureOrder = 2;
            lhs = LHSintegrator.create(s);
            K = lhs.compute();

        end
        function F = computeF(obj,displacementFun)

            s.type     = 'Elastic';
            s.scale    = 'MACRO';
            s.dim      = obj.getFunDims(displacementFun);
            s.BC       = obj.boundaryConditions;
            s.mesh     = obj.mesh;
            s.material = obj.material;
            s.globalConnec = obj.mesh.connec;
            RHSint = RHSintegrator.create(s);
            rhs = RHSint.compute();
            % NO Reactions where considered in this case
            F = rhs;
        end

        function dim = getFunDims(obj,displacementFun)
            d.ndimf  = displacementFun.ndimf;
            d.nnodes = size(displacementFun.fValues, 1);
            d.ndofs  = d.nnodes*d.ndimf;
            d.nnodeElem = obj.mesh.nnodeElem; % should come from interp..
            d.ndofsElem = d.nnodeElem*d.ndimf;
            dim = d;
        end

        function results = computeU(obj,K,F)
            s.type      = 'DIRECT';
            s.solverType = 'REDUCED'; %MIRAR COM INFLUEIX!!!
            s.solverMode = 'DISP';
            problemSolver = obj.createSolver(s);

            t.stiffness = K;
            t.forces = F;
            results = problemSolver.solve(t);

        end
        % FOR INTERNAL MATERIA CREATION
        % function createMaterial (obj)
        %     s.type = 'ISOTROPIC';
        %     j = Material.create()
        %
        % end

        function problemSolver = createSolver(obj,dataIn)
            sS.type      = dataIn.type;
            solver       = Solver.create(sS);
            s.solverType = dataIn.solverType;
            s.solverMode = dataIn.solverMode;
            s.solver     = solver;
            s.boundaryConditions = obj.boundaryConditions;
            s.BCApplier          = obj.createBCApplier();
            problemSolver    = ProblemSolver(s);
        end

        function BC = createBCApplier(obj)
            s.mesh = obj.mesh;
            s.boundaryConditions = obj.boundaryConditions;
            BC = BCApplier(s);
        end

    end

end