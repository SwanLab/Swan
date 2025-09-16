function LHS = IntegrateLHS(f,test,trial,mesh,quadOrder)
v = @(i) Test(test,i);
u = @(j) Test(trial,j);
f = @(i,j) f(u(j),v(i));
s.mesh  = mesh;
s.quadratureOrder = quadOrder;
lhs = LHSIntegrator(s);
LHS = lhs.compute(f,test,trial);
end