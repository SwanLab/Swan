classdef FilterPDE < handle
    
    properties (Access = private)
        mesh
        filteredField
        scale
        epsilon
        LHStype
        problemLHS
        LHS
        RHS
        bc
    end

    methods (Access = public)

        function obj = FilterPDE(cParams)
            obj.init(cParams);
            obj.computeBoundaryConditions();
            lhs     = obj.createProblemLHS(cParams);
            obj.LHS = decomposition(lhs);
        end

        function xF = compute(obj,fun,quadType)
            obj.computeRHS(fun,quadType);
            obj.solveFilter();
            xF = obj.filteredField;
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
            f                 = FilterPDEFactory(cParams);
            obj.LHStype       = f.LHStype;
            obj.scale         = f.scale;
            obj.mesh          = cParams.mesh;
            obj.filteredField = P1Function.create(obj.mesh,1); % trial will come from outside
            obj.epsilon       = cParams.mesh.computeMeanCellSize();
        end

        function computeRHS(obj,fun,quadType)
            int  = obj.computeRHSintegrator(quadType);
            test = obj.filteredField;
            rhs  = int.compute(fun,test);
            rhsR = obj.bc.fullToReducedVector(rhs);            
            obj.RHS = rhsR;
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

        function solveFilter(obj)
            s.type = 'DIRECT';
            solver = Solver.create(s);
            x   = solver.solve(obj.LHS,obj.RHS);
            xR  = obj.bc.reducedToFullVector(x);
            obj.filteredField.fValues = xR;  
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
            s.trial        = obj.filteredField;
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