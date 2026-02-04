




% Crear LagrangianFunction sobre una cohesiveMesh

cohesiveMesh = CohesiveMesh();
u   = LagrangianFunction.create(cohesiveMesh.mesh,2,'P1');

s.cohesiveMesh = cohesiveMesh;
s.u = u;
disp = CohesiveSeparationComputer(s);

