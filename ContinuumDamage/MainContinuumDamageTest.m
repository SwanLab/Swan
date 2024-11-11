clc;clear;close all

%load('TestForceTraction1Elem.mat')
%load('TestDisplacementTraction.mat')
cParams.mesh.meshLength = 1;
cParams.mesh.meshWidth = 1;
cParams.mesh.meshN = 1;
cParams.mesh.meshM = 1;



cParams.material.E = 1;
cParams.material.nu = 0;

cParams.bc.bcType = 'displacementBending'; %'FORCE'
cParams.bc.bcVal = 1;

cParams.solver.type = 'Elastic';
cParams.solver.solverType='REDUCED';
cParams.solver.solverMode = 'DISP';
cParams.solver.solverCase = 'DIRECT';
cParams.solver.scale = 'MACRO';

cParams.tol = 1e-10;

tester = TestingContinuumDamage(cParams);
data = tester.compute();
data.displacement.plot

tester.compareWithElasticProblem(data.displacement.fValues,uRef.fValues);