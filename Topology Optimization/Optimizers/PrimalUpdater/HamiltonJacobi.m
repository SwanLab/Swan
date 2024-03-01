classdef HamiltonJacobi < handle

    properties (Access = public)
        tau
    end

    properties (Access = private)
        phi
        filter
        epsilon
        velocity
    end

    methods (Access = public)
        function obj = HamiltonJacobi(cParams)
            obj.init(cParams);
            obj.setupFilter();
        end

        function x = update(obj,g,~)
            obj.computeVelocity(g);
            x = obj.computeNewLevelSet();
        end

        function computeFirstStepLength(obj,g,x,f)
            obj.tau = f*sqrt(norm(g)/norm(x));
        end

        function is = isTooSmall(obj)
            is = obj.tau < 1e-10;
        end

        function increaseStepLength(obj,f)
            obj.tau = f*obj.tau;
        end

        function decreaseStepLength(obj)
            obj.tau = obj.tau/2;
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.phi     = cParams.designVariable;
            obj.epsilon = cParams.designVariable.fun.mesh.computeMeanCellSize();
        end

        function computeVelocity(obj,g)
            s.mesh       = obj.phi.fun.mesh;
            s.fValues    = g;
            s.order      = 'P1';
            ss.fun       = LagrangianFunction(s);
            ss.uMesh     = obj.phi.getUnfittedMesh();
            unfFun       = UnfittedBoundaryFunction(ss);
            gFilter      = obj.filter.compute(unfFun,'QUADRATIC');
            V            = -gFilter.fValues;
            Vnorm        = max(abs(V(:)));
            obj.velocity = V/Vnorm;
        end

        function x = computeNewLevelSet(obj)
            x = obj.phi.fun.fValues;
            t = obj.tau;
            x = x - t*obj.velocity;
            x = obj.normalizeFunction(x);
        end

        function x = normalizeFunction(obj,x)
            m         = obj.phi.fun.mesh;
            s.fValues = x;
            s.mesh    = m;
            s.order   = 'P1';
            xFun      = LagrangianFunction(s);
            norm      = Norm.computeH1(m,xFun,obj.epsilon);
            xNorm     = sqrt(norm);
            x         = x/xNorm;
        end

        function setupFilter(obj)
            designVar           = obj.phi;
            s.mesh              = designVar.fun.mesh;
            s.designVarType     = designVar.type;
            s.scale             = 'MACRO';
            s.filterType        = 'PDE';
            s.quadType          = 'LINEAR';
            s.designVariable    = designVar;
            s.trial             = LagrangianFunction.create(s.mesh,1, 'P1');
            obj.filter          = Filter.create(s);
            obj.filter.updateEpsilon(obj.epsilon);
        end
    end

end