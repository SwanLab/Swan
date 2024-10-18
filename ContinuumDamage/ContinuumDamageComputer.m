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
        quadOrder

        ElasticFun
        ExternalWorkFun

        u
    end

    methods (Access = public)

        function obj = ContinuumDamageComputer(cParams)
            obj.init(cParams)
        end

        function results = compute(obj)
            bc = obj.boundaryConditions;
            uFun = LagrangianFunction.create(obj.mesh, obj.mesh.ndim, 'P1');
            uFun.fValues = obj.updateInitialDisplacement(bc,uFun);


            s.material = obj.material;
            s.u        = uFun;
            s.mesh     = obj.mesh;


            obj.ElasticFun = shFunc_Elastic(s);
            obj.ExternalWorkFun = shFunc_ExternalWork2(s);
            obj.u = uFun;

            K = obj.computeKInternal();

            Fext = obj.computeFExternal();
            Fint = obj.computeFInternal();
            F = Fext-Fint;

            results = obj.computeU(K,F);

        end
    end

    methods (Access = private)

        function init(obj,cParams)
            obj.quadOrder = 2;
            obj.mesh = cParams.mesh;
            obj.boundaryConditions = cParams.boundaryConditions;
            obj.material = cParams.material;
            obj.type      = cParams.type;
            obj.solverType = cParams.solverType;
            obj.solverMode = cParams.solverMode;
            obj.solverCase = cParams.solverCase;
            obj.scale = cParams.scale;
        end

        function K = computeKInternal(obj)%,dispFun)
            K = obj.ElasticFun.computeHessian(obj.quadOrder);
        end

        function F = computeFInternal(obj)%,displacementFun,K)
            Ftry = obj.ElasticFun.computeJacobian(obj.quadOrder);
            F = Ftry;
        end


        function F = computeFExternal(obj)
            fExt = obj.boundaryConditions.pointloadFun;
            Ftry = obj.ExternalWorkFun.computeGradient(obj.u,fExt,obj.quadOrder);
            F = Ftry;
        end

        function dim = getFunDims(obj,displacementFun)
            d.ndimf  = displacementFun.ndimf;
            d.nnodes = size(displacementFun.fValues, 1);
            d.ndofs  = d.nnodes*d.ndimf;
            d.nnodeElem = obj.mesh.nnodeElem; % should come from interp..
            d.ndofsElem = d.nnodeElem*d.ndimf;
            dim = d;
        end

        function u = computeU(obj,K,F)


            problemSolver = obj.createSolver();

            t.stiffness = K;
            t.forces = F;
            u = problemSolver.solve(t);

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


         function u = updateInitialDisplacement(obj,bc,uOld)
            restrictedDofs = bc.dirichlet_dofs;
            if isempty(restrictedDofs)
                u = uOld;
            else
                dirich = bc.dirichletFun;
                uVec = reshape(uOld.fValues',[uOld.nDofs 1]);
                dirichVec = reshape(dirich.fValues',[dirich.nDofs 1]);

                uVec(restrictedDofs) = dirichVec(restrictedDofs);
                u = reshape(uVec,[flip(size(uOld.fValues))])';
            end
         end
    end
end