classdef FilterPDE < handle
    
    properties (Access = private)
        mesh
        trial
        epsilon
        LHStype
    end

    properties (Access = private)
        problemLHS
        LHS
        RHS
        % bc
    end

    methods (Access = public)
        function obj = FilterPDE(cParams)
            obj.init(cParams);
            % obj.computeBoundaryConditions(cParams);
            obj.createProblemLHS(cParams);
            obj.computeLHS();
        end

        function xF = compute(obj,fun,quadType)
            xF = LagrangianFunction.create(obj.mesh, fun.ndimf, obj.trial.order);
            obj.computeRHS(fun,quadType);
            obj.solveFilter();
            xF.setFValues(obj.trial.fValues);
        end

        function obj = updateEpsilon(obj,epsilon)
            if obj.hasEpsilonChanged(epsilon)
                obj.epsilon = epsilon;
                obj.computeLHS();
            end
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.trial   = LagrangianFunction.create(cParams.mesh, cParams.trial.ndimf, cParams.trial.order);
            obj.LHStype = cParams.LHStype;
            obj.mesh    = cParams.mesh;
            obj.epsilon = cParams.mesh.computeMeanCellSize();
        end

        function computeBoundaryConditions(obj,cParams)
            s                 = cParams;
            s.bc{1}.dirichlet = [];
            s.bc{1}.pointload = [];
            s.bc{1}.ndimf     = 1;
            s.bc{1}.ndofs     = [];
            s.ndofs           = obj.mesh.nnodes;
            obj.bc            = BoundaryConditions(s);
            obj.bc.compute();
        end

        function createProblemLHS(obj,s)
            s.trial        = obj.trial;
            s.mesh         = obj.mesh;
            s.type         = obj.LHStype;
            obj.problemLHS = LHSintegrator.create(s);
        end

        function computeLHS(obj)
            lhs     = obj.problemLHS.compute(obj.epsilon);
            % lhs     = obj.bc.fullToReducedMatrix(lhs); % its the same
            obj.LHS = decomposition(lhs);
        end

        function computeRHS(obj,fun,quadType)
            switch class(fun)
                case {'UnfittedFunction','UnfittedBoundaryFunction'}
                    s.mesh = fun.unfittedMesh;
                    s.type = 'Unfitted';
                otherwise
                    s.mesh = obj.mesh;
                    s.type = 'ShapeFunction';
            end
            s.quadType = quadType;
            int        = RHSintegrator.create(s);
            test       = obj.trial;
            rhs        = int.compute(fun,test);
            % rhsR       = obj.bc.fullToReducedVector(rhs);
            rhsR       = rhs; % its the same
            obj.RHS    = rhsR;
        end

        function solveFilter(obj)
            s.type = 'DIRECT';
            solver = Solver.create(s);
            x      = solver.solve(obj.LHS,obj.RHS);
            % xR     = obj.bc.reducedToFullVector(x);
            xR = x; % its the same
            obj.trial.setFValues(reshape(xR',obj.trial.ndimf,[])');
        end

        function itHas = hasEpsilonChanged(obj,eps)
            var   = abs(eps - obj.epsilon)/eps;
            itHas = var > 1e-15;
        end
    end
end