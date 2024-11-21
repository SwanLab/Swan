clc;clear;close all

%load('TestForceTraction1Elem.mat')
%load('TestDisplacementTraction.mat')
cParams.mesh.meshLength = 1;
cParams.mesh.meshWidth = 1;
cParams.mesh.meshN = 10;
cParams.mesh.meshM = 50;



cParams.material.E = 3e4;
cParams.material.nu = 0.3;

cParams.bc.bcType = 'forceTraction'; %'FORCE'
cParams.bc.bcVal = 1;

cParams.solver.type = 'Elastic';
cParams.solver.solverType='REDUCED';
cParams.solver.solverMode = 'DISP';
cParams.solver.solverCase = 'DIRECT';
cParams.solver.scale = 'MACRO';

cParams.tol = 1e-10;

tester = TestingContinuumDamage(cParams);
data = tester.compute();


data.displacement.plot()
data.reactions.plot()
data.damage.plot(data.displacement.mesh)


tester.compareWithElasticProblem(data.displacement.fValues,uRef.fValues);