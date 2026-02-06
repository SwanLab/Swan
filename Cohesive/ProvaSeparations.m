




% Crear LagrangianFunction sobre una cohesiveMesh

cohesiveMesh = CohesiveMesh();
u   = LagrangianFunction.create(cohesiveMesh.mesh,2,'P1');
% fixar valors de u a cohesius per comprovar

s.cohesiveMesh = cohesiveMesh;
s.u = u;
s.ndimf = 2;
disp = CohesiveSeparationComputer(s);


disp.compute(u);

