
clear;
clc;
close all;

% Crear LagrangianFunction sobre una cohesiveMesh
mesh = UnitQuadMesh(2,2);

ymin = min(mesh.coord(:,2));

s.baseMesh    = mesh;
s.isFractured = @(coord) abs(coord(:,2) - ymin) <=1e-10;
cohesiveMesh = CohesiveMesh(s);


u  = LagrangianFunction.create(cohesiveMesh.mesh,2,'P1');
%% Comprovacions de separacions
fValues = u.fValues;

comprovacio = 1;

if comprovacio == 1
    separacio = -0.1;
    fValues(1,2) = separacio;
    fValues(4,2) = separacio;
    fValues(7,2) = separacio;

    u.setFValues(fValues);

elseif comprovacio == 2
    fValues(1,2) = -0.1;
    fValues(5,2) = -0.075;
    fValues(9,2) = -0.05;
    fValues(13,2) = -0.025;
    u.setFValues(fValues);

elseif comprovacio == 3
    fValues(1,2) = -0.025;
    fValues(5,2) = -0.05;
    fValues(9,2) = -0.075;
    fValues(13,2) = -0.1;
    u.setFValues(fValues);

elseif comprovacio == 4
    separacio = 0.1;
    fValues(1,2) = separacio;
    fValues(5,2) = separacio;
    fValues(9,2) = separacio;
    fValues(13,2) = separacio;
    u.setFValues(fValues);
end


%%

s.cohesiveMesh = cohesiveMesh;
s.uFun = u;
s.ndimf = 2;
separator = Jump(s);


jump = separator.jumpFun;
xV = [-1 1];

separator.computeShapeFunctions(xV);












        s.lawType               = 'Cubic';
        s.normalCharLength      = 0.01;
        s.tangencialCharLength  = 0.01;
        s.sigmaMax              = 1;

constitutiveLaw = CohesiveTractionSeparation(s);
constitutiveLaw.evaluate(xV,jump)