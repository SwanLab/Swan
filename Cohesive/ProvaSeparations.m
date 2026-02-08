
clear;
clc;
close all;

% Crear LagrangianFunction sobre una cohesiveMesh
s.baseMesh = UnitQuadMesh(4,4);
cohesiveMesh = CohesiveMesh(s);
u   = LagrangianFunction.create(cohesiveMesh.mesh,2,'P1');


%% Comprovacions de separacions
fValues = u.fValues;

comprovacio = 4;

if comprovacio == 1
    separacio = -0.1;
    fValues(1,2) = separacio;
    fValues(5,2) = separacio;
    fValues(9,2) = separacio;
    fValues(13,2) = separacio;
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
    s.u = u;
    s.ndimf = 2;
separator = CohesiveSeparationComputer(s);

separator.compute(u);

jump = separator.fun;

disp('Tangencial - Normal')
disp(jump.fValues);








