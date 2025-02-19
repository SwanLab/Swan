mesh = UnitQuadMesh(10,10);

u = LagrangianFunction.create(mesh,2,'P2');
v = LagrangianFunction.create(mesh,2,'P2');

p = LagrangianFunction.create(mesh,2,'P1');
q = LagrangianFunction.create(mesh,2,'P1');


grad_u = Grad(u);
grad_v = Grad(v);
div_u = Divergence(u);

A = Integrator.compute(DP(grad_u, grad_v), mesh, 2);
a = Voigt(DP(DP(u, grad_u),v));
N = Integrator.compute(a, mesh, 2);
B = Integrator.compute(DP(q,div_u), mesh, 2);