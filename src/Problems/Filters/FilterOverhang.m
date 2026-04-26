classdef FilterOverhang < handle
    
    properties (Access = private)
        mesh
        trial
        k
        theta
        epsilon
    end

    properties (Access = private)
        chiN
        chiNOld
    end

    methods (Access = public)
        function obj = FilterOverhang(cParams)
            obj.init(cParams);
        end

        function xF = compute(obj,fun,q)
            xF = copy(obj.trial);
            obj.chiN = obj.createRHSShapeFunction(fun,q);
            if isempty(obj.chiNOld) || norm(obj.chiNOld-obj.chiN)/norm(obj.chiN)>=1e-6
                obj.chiNOld      = obj.chiN;
                s.designVariable = obj.createDesignVariable(xF);
                s.cost           = obj.createCost(fun);
                s.monitoring     = false;
                s.lb             = 0;
                s.ub             = 1;
                s.maxIter        = 1000;
                opt              = OptimizerProjectedGradient(s);
                opt.solveProblem();
            end
        end

        function obj = updateEpsilon(obj,epsilon)
            if obj.hasEpsilonChanged(epsilon)
                obj.epsilon = epsilon;
            end
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.trial   = LagrangianFunction.create(cParams.mesh, cParams.trial.ndimf, cParams.trial.order);
            obj.mesh    = cParams.mesh;
            obj.k       = cParams.senseVector;
            obj.theta   = deg2rad(cParams.ovAngleDeg);
            obj.epsilon = cParams.mesh.computeMeanCellSize();
        end

        function RHS = createRHSShapeFunction(obj,fun,quadType)
            f   = @(v) DP(fun,v);
            RHS = IntegrateRHS(f,obj.trial,obj.trial.mesh,'Domain',quadType);
        end

        function dens = createDesignVariable(obj,xF)
            s.fun      = xF;
            s.mesh     = obj.mesh;
            s.type     = 'Density';
            s.plotting = false;
            dens       = DesignVariable.create(s);
        end

        function c = createCost(obj,fun)
            s.shapeFunctions{1} = obj.createPotential(fun);
            s.weights           = 1;
            s.Msmooth           = obj.createMassMatrix();
            c                   = Cost(s);
        end

        function p = createPotential(obj,fun)
            s.mesh    = obj.mesh;
            s.chi     = fun;
            s.epsilon = obj.epsilon;
            s.k       = obj.k;
            s.theta   = obj.theta;
            s.trial   = obj.trial;
            p         = OverhangPotential(s);
        end

        function M = createMassMatrix(obj)
            n = obj.mesh.nnodes;
            h = obj.mesh.computeMinCellSize();
            M = h^2*sparse(1:n,1:n,ones(1,n),n,n);
        end

        function itHas = hasEpsilonChanged(obj,eps)
            var   = abs(eps - obj.epsilon)/eps;
            itHas = var > 1e-15;
        end
    end
end