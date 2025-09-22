function RHS = IntegrateRHS(f,test,mesh,quadOrder)
if nargin < 4 || isempty(quadOrder)
    quadOrder = 2;
end
s.test  = test;
s.mesh  = mesh;
s.quadratureOrder = quadOrder;
rhs = RHSIntegrator(s);
RHS = rhs.compute(f);
end