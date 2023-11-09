classdef HamiltonJacobi < handle

    properties (Access = public)
        tau
    end

    properties (Access = private)
        phi
        filter
        scalar_product
        velocity
    end

    methods (Access = public)
        
        function obj = HamiltonJacobi(cParams)
            obj.init(cParams);
            obj.setupFilter(obj.scalar_product.epsilon,obj.phi);
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
            obj.phi            = cParams.designVar;
            obj.scalar_product = cParams.uncOptimizerSettings.scalarProductSettings;
        end

        function computeVelocity(obj,g)
            s.mesh       = obj.phi.mesh;
            s.fValues    = g;
            ss.fun       = P1Function(s);
            ss.uMesh     = obj.phi.getUnfittedMesh();
            unfFun       = UnfittedBoundaryFunction(ss);
            V            = -obj.filter.compute(unfFun,'QUADRATICMASS');
            Vnorm        = max(abs(V(:)));
            obj.velocity = V/Vnorm;
        end

        function x = computeNewLevelSet(obj)
            x = obj.phi.value;
            t = obj.tau;
            x = x - t*obj.velocity;
            x = obj.normalizeFunction(x);
        end

        function x = normalizeFunction(obj,x)
            norm2 = obj.scalar_product.computeSP(x,x);
            xNorm = sqrt(norm2);
            x = x/xNorm;
        end

        function setupFilter(obj,e,designVar)
            set                 = SettingsFilter('paramsFilter_PDE_Boundary.json');
            s                   = set.femSettings;
            s.mesh              = designVar.mesh;
            s.designVarType     = designVar.type;
            s.scale             = 'MACRO';
            s.filterType        = set.filterType;
            s.quadType          = 'LINEAR';
            s.designVariable    = designVar;
            obj.filter          = Filter.create(s);
            obj.filter.updateEpsilon(e);
        end


    end


end