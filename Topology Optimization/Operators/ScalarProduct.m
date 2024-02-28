classdef ScalarProduct < handle

    methods (Access = public, Static)

        function sp = computeL2(m,f,g)
            int = Integrator.create('ScalarProduct',m,'QUADRATIC');
            sp = int.compute(f,g);
        end

        function sp = computeH1(m,f,g,eps)
            quadOrder = 'QUADRATIC';
            q = Quadrature.set(m.type);
            q.computeQuadrature(quadOrder);
            Df  = Grad(f);
            Dg  = Grad(g);
            int = Integrator.create('ScalarProduct',m,quadOrder);
            spM = int.compute(f,g);
            spK = int.compute(Df,Dg);
            sp  = spM + eps^2*spK;
        end

    end

end
