clc;clear;close all

%load('TestForceTraction1Elem.mat')
%load('TestDisplacementTraction.mat')
%cParams.mesh.name = 'CD_Mesh';

% type = '1Dtrac';
% 
% switch type
%     case 'SEMtrac'
%         cParams.mesh.name = 'PF_SENtraction0_0025';
%         cParams.bc.bcType = 'SEMtraction'; 
%     case 'SEMmix'
%         cParams.mesh.name = 'PF_SENmixed0_0025';
%         cParams.bc.bcType = 'SEMmixed'; 
%     case 'SEMshear'
%         cParams.mesh.name = 'PF_SENshear0_0025';
%         cParams.bc.bcType = 'SEMshear'; 
%     case '1Dtrac'
%         cParams.mesh.meshLength = 1;
%         cParams.mesh.meshWidth = 1;
%         cParams.mesh.meshN = 10;
%         cParams.mesh.meshM = 10;
%         cParams.bc.bcType = 'displacementTraction';
% end



cParams.mesh.meshLength = 1;
cParams.mesh.meshWidth = 1;
cParams.mesh.meshN = 10;
cParams.mesh.meshM = 10;

cParams.bc.bcType = 'displacementTraction'; %'FORCE'
cParams.bc.bcValueSetLoading = 1e-10:1e-3:0.2;
cParams.bc.bcValueSetUnLoading = [];

cParams.material.E = 210;
cParams.material.nu = 0.3;

cParams.solver.type = 'Elastic';
cParams.solver.solverType='REDUCED';
cParams.solver.solverMode = 'DISP';
cParams.solver.solverCase = 'DIRECT';
cParams.solver.scale = 'MACRO';

cParams.tol = 1e-8;
cParams.H = 0.5;
cParams.r0 = 1/sqrt(6);
cParams.r1 = 2;

tester = TestingContinuumDamage(cParams);
data = tester.compute();

data.displacement.field.plot
data.damage.field.plot
colorbar
caxis([0 1-cParams.H])

figure()
plot(data.displacement.value,data.damage.maxValue)
title('Damage-Displacement')

figure()
plot(data.displacement.value,data.reaction)
title('Force-Displacement')

figure()
plot(data.displacement.value,data.totalEnergy)
title('Energy - Displacement')

figure()
plot(data.displacement.value,data.damagedMaterial)
title('Material - Displacement')

%tester.compareWithElasticProblem(data.displacement.fValues,uRef.fValues);