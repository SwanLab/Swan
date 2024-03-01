classdef ScalarProduct < handle

    methods (Access = public, Static)

        function sp = computeL2(m,f,g)
            fg = DDP(f,g);
            sp = Integrator.compute(fg,m,'QUADRATIC');
        end

        function sp = computeH1(m,f,g,eps)
            quadOrder = 'QUADRATIC';
            q = Quadrature.set(m.type);
            q.computeQuadrature(quadOrder);
            Df   = Grad(f);
            Dg   = Grad(g);
            fg   = DDP(f,g);
            DfDg = DDP(Df,Dg);
            spM = Integrator.compute(fg,m,quadOrder);
            spK = Integrator.compute(DfDg,m,quadOrder);
            sp  = spM + eps^2*spK;
        end

    end

end
