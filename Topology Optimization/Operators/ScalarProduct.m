classdef ScalarProduct < handle

    methods (Access = public, Static)

        function sp = computeL2(m,f,g)
            %fg = DDP(f,g);
            fg = f.*g;
            sp = Integrator.compute(fg,m,2);
        end

        function sp = computeH1(m,f,g,eps)
            quadOrder = 2;
            spM = ScalarProduct.computeL2(m,f,g);           
            Df   = Grad(f);
            Dg   = Grad(g);
            DfDg = DP(Df,Dg);
            spK = Integrator.compute(DfDg,m,quadOrder);
            sp  = spM + eps^2*spK;
        end

    end

end
