clc;clear;close all

%load('TestForceTraction1Elem.mat')
%load('TestDisplacementTraction.mat')
cParams.mesh.meshLength = 1;
cParams.mesh.meshWidth = 1;
cParams.mesh.meshN = 30;
cParams.mesh.meshM = 30;



cParams.material.E = 210;
cParams.material.nu = 0.3;

cParams.bc.bcType = 'displacementTraction'; %'FORCE'
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



TotalReac = sum(data.reactions.fValues(:,2));
data.damage.plot(data.displacement.mesh)
damageFun = data.damage.project('P1D',data.displacement.mesh);


tester.compareWithElasticProblem(data.displacement.fValues,uRef.fValues);