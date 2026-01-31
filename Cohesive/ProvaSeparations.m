




% Crear LagrangianFunction sobre una cohesiveMesh

cohesiveMesh = CohesiveMesh();
u   = LagrangianFunction.create(cohesiveMesh.mesh,2,'P1');


disp = CohesiveSeparationComputer(u,cohesiveMesh);