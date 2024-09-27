classdef ContinuumDamageComputer < handle

    properties (Access = private)
        mesh
        boundaryConditions
        material

        type
        scale
        solverType
        solverMode
        solverCase
    end

    methods (Access = public)

        function obj = ContinuumDamageComputer(cParams)
            obj.init(cParams)
        end

        function results = compute(obj)
            displacementFun = LagrangianFunction.create(obj.mesh, obj.mesh.ndim, 'P1');
            K = obj.computeK(displacementFun);
            F = obj.computeF(displacementFun,K);
            results = obj.computeU(K,F);
        end
    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh = cParams.mesh;
            obj.boundaryConditions = cParams.boundaryConditions;
            obj.material = cParams.material;
            obj.type      = cParams.type;
            obj.solverType = cParams.solverType; %MIRAR COM INFLUEIX!!!
            obj.solverMode = cParams.solverMode;
            obj.solverCase = cParams.solverCase;
            obj.scale = cParams.scale;
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
        function F = computeF(obj,displacementFun,K)

            s.type     = obj.type;
            s.scale    = obj.scale;
            s.dim      = obj.getFunDims(displacementFun);
            s.BC       = obj.boundaryConditions;
            s.mesh     = obj.mesh;
            s.material = obj.material;
            s.globalConnec = obj.mesh.connec;
            RHSint = RHSintegrator.create(s);
            rhs = RHSint.compute();
          
            if strcmp(obj.solverType,'REDUCED')
                R = RHSint.computeReactions(K);
                F = rhs+R;
            else
                F = rhs;
            end

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

            s.solverType = obj.solverType; %MIRAR COM INFLUEIX!!!
            s.solverMode = obj.solverMode;
            problemSolver = obj.createSolver();

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

        function problemSolver = createSolver(obj)
            sS.type      = obj.solverCase;
            solver       = Solver.create(sS);
            s.solverType = obj.solverType;
            s.solverMode = obj.solverMode;
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