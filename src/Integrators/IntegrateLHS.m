function LHS = IntegrateLHS(f,test,trial,mesh,quadOrder)
s.trial = trial;
s.test  = test;
s.mesh  = mesh;
s.quadratureOrder = quadOrder;
lhs = LHSIntegrator(s);
LHS = lhs.compute(f);
end