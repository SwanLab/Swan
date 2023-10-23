classdef FilterPDE < handle
    
    properties (Access = private)
        mesh
        trial
        scale
        epsilon
        LHStype
        problemLHS
        LHS
        bc
    end

    methods (Access = public)

        function obj = FilterPDE(cParams)
            obj.init(cParams);
            obj.computeBoundaryConditions();
            lhs     = obj.createProblemLHS(cParams);
            obj.LHS = decomposition(lhs);
        end

        function xReg = compute(obj,fun,quadType)
            RHS               = obj.computeRHS(fun,quadType);
            xR                = obj.solveFilter(RHS);
            obj.trial.fValues = xR;
            xReg              = obj.trial;
        end

        function obj = updateEpsilon(obj,epsilon)
            if obj.hasEpsilonChanged(epsilon)
                obj.epsilon = epsilon;
                lhs         = obj.computeLHS();
                obj.LHS     = decomposition(lhs);
            end
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh  = cParams.mesh;
            obj.trial = P1Function.create(obj.mesh,1);
            obj.scale = cParams.scale;
            if isfield(cParams,'LHStype')
                obj.LHStype = cParams.LHStype;
            else
                obj.LHStype = 'DiffReactNeumann';
            end
            obj.epsilon = cParams.mesh.computeMeanCellSize();
        end

        function RHS = computeRHS(obj,fun,quadType)
            int  = obj.computeRHSintegrator(quadType);
            test = P1Function.create(obj.mesh, 1);
            RHS  = int.compute(fun,test);
        end

        function int = computeRHSintegrator(obj,quadType)
            s.mesh     = obj.mesh;
            s.type     = 'ShapeFunction';
            s.quadType = quadType;
            int        = RHSintegrator.create(s);
        end

        function itHas = hasEpsilonChanged(obj,eps)
            if isempty(obj.epsilon)
                obj.epsilon = 0;
            end
            var = abs(eps - obj.epsilon)/eps;
            itHas = var > 1e-15;
        end

        function x_reg = solveFilter(obj,RHS)
            RHS    = obj.bc.fullToReducedVector(RHS);
            s.type = 'DIRECT';
            Solv   = Solver.create(s);
            x      = Solv.solve(obj.LHS,RHS);
            x_reg  = obj.bc.reducedToFullVector(x);
        end

        function computeBoundaryConditions(obj)
            s.scale           = obj.scale;
            s.mesh            = obj.mesh;
            s.bc{1}.dirichlet = [];
            s.bc{1}.pointload = [];
            s.bc{1}.ndimf     = 1; % periodic BCs
            s.bc{1}.ndofs     = [];
            s.ndofs           = obj.mesh.nnodes;
            obj.bc            = BoundaryConditions(s);
            obj.bc.compute();
        end

        function lhs = createProblemLHS(obj,s)
            s.trial        = obj.trial;
            s.mesh         = obj.mesh;
            s.type         = obj.LHStype;
            obj.problemLHS = LHSintegrator.create(s);
            lhs            = obj.computeLHS();
        end

        function lhs = computeLHS(obj)
            lhs = obj.problemLHS.compute(obj.epsilon);
            lhs = obj.bc.fullToReducedMatrix(lhs);
        end

    end

end