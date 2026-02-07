




% Crear LagrangianFunction sobre una cohesiveMesh

cohesiveMesh = CohesiveMesh();
u   = LagrangianFunction.create(cohesiveMesh.mesh,2,'P1');


%% Comprovacions de separacions


fValues = u.fValues;


disp = -0.1;

fValues(1,2) = disp;
fValues(5,2) = disp;
fValues(9,2) = disp;
fValues(13,2) = disp;

u.setFValues(fValues);








s.cohesiveMesh = cohesiveMesh;
s.u = u;
s.ndimf = 2;
disp = CohesiveSeparationComputer(s);


disp.compute(u);










