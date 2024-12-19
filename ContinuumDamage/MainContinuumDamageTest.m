clc;clear;close all

%load('TestForceTraction1Elem.mat')
%load('TestDisplacementTraction.mat')
cParams.mesh.meshLength = 1;
cParams.mesh.meshWidth = 1;
cParams.mesh.meshN = 10;
cParams.mesh.meshM = 10;

%cParams.mesh.name = 'CD_Mesh';

cParams.material.E = 210;
cParams.material.nu = 0.3;

cParams.bc.bcType = 'displacementTraction'; %'FORCE'
cParams.bc.bcValueSet = [0.02,0.021];


cParams.solver.type = 'Elastic';
cParams.solver.solverType='REDUCED';

cParams.solver.solverMode = 'DISP';
cParams.solver.solverCase = 'DIRECT';
cParams.solver.scale = 'MACRO';

cParams.tol = 1e-10;

tester = TestingContinuumDamage(cParams);
data = tester.compute();

data.displacement.field.plot
data.damage.field.plot

figure()
plot(data.displacement.value,data.damage.maxValue)
title('Damage-Displacement')

figure()
plot(data.displacement.value,data.reaction)
title('Force-Displacement')

figure()
plot(data.displacement.value,data.totalEnergy)
title('Energy - Displacement')

tester.compareWithElasticProblem(data.displacement.fValues,uRef.fValues);