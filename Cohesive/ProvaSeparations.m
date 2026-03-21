
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
    fValues(1,2) = separacio *0.25;
    fValues(5,2) = separacio * 0.5;
    fValues(9,2) = separacio * 1;
    fValues(13,2) = separacio * 5;
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
g.jumpCrit  = 0.02;
g.K = 1e15;
tractionSeparation = CohesiveTractionSeparation(g);

% Calcular tractions
traction = tractionSeparation.computeFunction(jump);
derivative = tractionSeparation.computeDerivative(jump);
derivative.evaluate([-1,1]);

% Crear material 
k.type    = 'ISOTROPIC';
k.ptype   = 'ELASTIC';
k.ndim    = cohesiveMesh.mesh.ndim;
k.young   = ConstantFunction.create(100,cohesiveMesh.mesh);
k.poisson = ConstantFunction.create(0.33, cohesiveMesh.mesh);
k.cohesiveMesh = cohesiveMesh;
material    = Material.create(k);

% Funcionals
p.tractionSeparation = tractionSeparation;
p.material     = material;
p.cohesiveMesh = cohesiveMesh;
p.test = LagrangianFunction.create(cohesiveMesh.mesh,2,'P1');
funcional = CohesiveFunctional(s);


