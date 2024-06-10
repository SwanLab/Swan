classdef ScalarProduct < handle

    methods (Access = public, Static)

        function sp = computeL2(m,f,g)
            fg = f.*g;
            sp = Integrator.compute(fg,m,2);
        end

        function sp = computeH1(m,f,g,eps)
            quadOrder = 2;
            q = Quadrature.create(m,quadOrder);
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
