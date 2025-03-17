function sp = ScalarProduct(f,g,type,eps)
switch type
    case 'L2'
        sp = computeL2(f,g);
    case 'H1'
        sp = computeH1(f,g,eps);
end
end

function sp = computeL2(f,g)
fg = f.*g;
sp = Integrator.compute(fg,fg.mesh,2);
end

function sp = computeH1(f,g,eps)
quadOrder = 2;
spM  = computeL2(f,g);
Df   = Grad(f);
Dg   = Grad(g);
DfDg = DP(Df,Dg);
spK = Integrator.compute(DfDg,DfDg.mesh,quadOrder);
sp  = spM + eps^2*spK;
end

