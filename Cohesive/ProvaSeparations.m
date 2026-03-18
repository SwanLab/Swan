
clear;
clc;
close all;

% Crear LagrangianFunction sobre una cohesiveMesh
mesh = UnitQuadMesh(3,3);

% Condicions
ymin = min(mesh.coord(:,2));

% Crear Cohesive MEsh
s.baseMesh    = mesh;
s.isFractured = @(coord) abs(coord(:,2) - ymin) <=1e-10;
cohesiveMesh = CohesiveMesh(s);

% Inicialitzar functio de desplaçaments
u  = LagrangianFunction.create(cohesiveMesh.mesh,2,'P1');

% Comprovacions de separacions
fValues = u.fValues;
    separacio = -0.1;
    fValues(1,2) = separacio;
    fValues(5,2) = separacio;
    fValues(9,2) = separacio;
    fValues(13,2) = separacio;
    u.setFValues(fValues);


% Inicialització classe Jump
f.cohesiveMesh = cohesiveMesh;
f.uFun = u;
f.ndimf = 2;
jump = Jump(f);

% Comprovació Funcions de Forma
xV = [-1 1];
Bc = jump.computeShapeFunctions(xV);


% ================================================== %

% Inicialitzar Traction Separation

g.lawType = 'TractionBiliniarUncoupled';
g.jumpFinal = 0.1;
g.jumpCrit  = 0.05;
g.K = 1000;
tractionSeparation = CohesiveTractionSeparation(g);

tractionSeparation.computeFunction(jump.jumpFun);