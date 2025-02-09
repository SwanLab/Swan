classdef H1ScalarProduct < handle

    properties (Access = private)
        integrator
        quadrature
    end

    methods (Access = public)
        function obj = H1ScalarProduct(mesh)
            obj.createIntegrator(mesh);
            obj.createQuadrature(mesh);
        end

        function sp = compute(obj,f,g,epsilon)
            Df  = f.evaluateGradient(obj.quadrature.posgp);
            Dg  = g.evaluateGradient(obj.quadrature.posgp);
            spM = obj.integrator.compute(f,g);
            spK = obj.integrator.compute(Df,Dg);
            sp  = spM+epsilon^2*spK;
        end
    end

    methods (Access = private)
        function createIntegrator(obj,m)
            q.mesh         = m;
            q.quadType     = 'QUADRATIC';
            q.type         = 'ScalarProduct';
            int            = Integrator.create(q);
            obj.integrator = int;
        end

        function createQuadrature(obj,m)
            q = Quadrature.set(m.type);
            q.computeQuadrature('QUADRATIC');
            obj.quadrature = q;
        end
    end
end