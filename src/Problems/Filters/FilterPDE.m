classdef FilterPDE < handle
    
    properties (Access = private)
        mesh
        trial
        epsilon
        LHSint
    end

    properties (Access = private)
        LHS
        RHS
        bc
    end

    methods (Access = public)
        function obj = FilterPDE(cParams)
            obj.init(cParams);
            obj.computeBoundaryConditions(cParams);
            obj.computeLHS();
        end

        function xF = compute(obj,fun,quadType)
            xF = LagrangianFunction.create(obj.mesh, fun.ndimf, obj.trial.order);
            obj.computeRHS(fun,quadType);
            obj.solveFilter();
            xF.setFValues(full(obj.trial.fValues));
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
            obj.LHSint  = cParams.LHSint;
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
            obj.bc            = BoundaryConditionsStokes(s); % This does not make sense regarding clean code. Once ProblemSolver is refactored by Raul and Alex, this will be changed. Speak with Jose when you read this.
            obj.bc.compute();
        end


        function computeLHS(obj)
            e   = obj.epsilon;
            vF  = obj.trial;
            uF  =  obj.trial;
            lhs = IntegrateLHS(@(u,v) obj.LHSint.domain(e,u,v),vF,uF,obj.mesh,'Domain');
            if ~isempty(obj.LHSint.boundary)
                lhs = lhs + IntegrateLHS(@(u,v) obj.LHSint.boundary(e,u,v),vF,uF,obj.mesh,'Boundary');
            end
            lhs     = obj.bc.fullToReducedMatrix(lhs);
            obj.LHS = decomposition(lhs);
        end

        function computeRHS(obj,fun,quadType)
            switch class(fun)
                case {'UnfittedFunction','UnfittedBoundaryFunction'}
                    s.mesh = fun.unfittedMesh;
                    s.type = 'Unfitted';
                    s.quadType = quadType;
                    int        = RHSIntegrator.create(s);
                    obj.RHS    = int.compute(fun,obj.trial);
                otherwise
                    f = @(v) DP(v,fun);
                    obj.RHS = IntegrateRHS(f,obj.trial,obj.trial.mesh,quadType);   
            end            
             obj.RHS     = obj.bc.fullToReducedVector(obj.RHS);
        end

        function solveFilter(obj)
            s.type = 'DIRECT';
            solver = Solver.create(s);
            x      = solver.solve(obj.LHS,obj.RHS);
            xR     = obj.bc.reducedToFullVector(x);
            obj.trial.setFValues(reshape(xR',obj.trial.ndimf,[])');
        end

        function itHas = hasEpsilonChanged(obj,eps)
            var   = abs(eps - obj.epsilon)/eps;
            itHas = var > 1e-15;
        end
    end
end