clear;clc;
cParams.meshLength = 1;
cParams.meshWidth = 1;
cParams.meshN = 1;
cParams.meshM = 1;
cParams.E = 210;
cParams.nu = 0.3;
cParams.bcVal = 1;
cParams.ndim = 2;
tol = 1e-15;
type = 'DISP';

cParams.type = 'Elastic';
cParams.solverType='REDUCED';
cParams.solverMode = 'DISP';
cParams.solverCase = 'DIRECT';
cParams.scale = 'MACRO';

main = TestingContinuumDamage(cParams,tol,type);