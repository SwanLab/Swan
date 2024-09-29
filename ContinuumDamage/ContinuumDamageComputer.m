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
    end

    methods (Access = public)

        function obj = ContinuumDamageComputer(cParams)
            obj.init(cParams)
        end

        function results = compute(obj)
            displacementFun = LagrangianFunction.create(obj.mesh, obj.mesh.ndim, obj.quadOrder);
            s.material = obj.material;
            
            s.u        = displacementFun;
            s.mesh     = obj.mesh;
            
            obj.ElasticFun = shFunc_Elastic(s);
            %aa = obj.ElasticFun.computeFunction(1)

            K = obj.computeK();
            F = obj.computeF();
            
            
            results = obj.computeU(K,F);

        end
    end

    methods (Access = private)

        function init(obj,cParams)
            obj.quadOrder = 'P1';
            obj.mesh = cParams.mesh;
            obj.boundaryConditions = cParams.boundaryConditions;
            obj.material = cParams.material;
            obj.type      = cParams.type;
            obj.solverType = cParams.solverType;
            obj.solverMode = cParams.solverMode;
            obj.solverCase = cParams.solverCase;
            obj.scale = cParams.scale;
        end

        function K = computeK(obj)%,dispFun)

            % ndimf = dispFun.ndimf;
            % s.type     = 'ElasticStiffnessMatrix';
            % s.mesh     = obj.mesh;
            % s.test     = LagrangianFunction.create(obj.mesh,ndimf, 'P1');
            % s.trial    = dispFun;
            % s.material = obj.material;
            % s.quadratureOrder = 2;
            % lhs = LHSintegrator.create(s);
            % K = lhs.compute();
            ord = obj.convertOrder();
            K = obj.ElasticFun.computeJacobian(ord);
        end
        function F = computeF(obj)%,displacementFun,K)

            % s.type     = obj.type;
            % s.scale    = obj.scale;
            % s.dim      = obj.getFunDims(displacementFun);
            % s.BC       = obj.boundaryConditions;
            % s.mesh     = obj.mesh;
            % s.material = obj.material;
            % s.globalConnec = obj.mesh.connec;
            % RHSint = RHSintegrator.create(s);
            % rhs = RHSint.compute();
            %
            % if strcmp(obj.solverType,'REDUCED')
            %     R = RHSint.computeReactions(K);
            %     F = rhs+R;
            % else
            %     F = rhs;
            % end
            ord = obj.convertOrder();
            F = obj.ElasticFun.computeHessian(ord);
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

            s.solverType = obj.solverType; %MIRAR COM INFLUEIX!!!
            s.solverMode = obj.solverMode;
            problemSolver = obj.createSolver(s);

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

        function ord = convertOrder(obj,order)
            if ischar(obj.quadOrder)
                switch obj.quadOrder
                    case 'P0'
                        ord = 0;
                    case 'P1'
                        ord = 1;
                    case 'P2'
                        ord = 2;
                    case 'P3'
                        ord = 3;
                end
            else
                ord = order;
            end
        end

    end
end