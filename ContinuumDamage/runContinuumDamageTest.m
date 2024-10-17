clc;clear;

%load('ContinuumDamageTestInput.mat')

cParams.meshLength = 1;
cParams.meshWidth = 1;
cParams.meshN = 2;
cParams.meshM = 2;
cParams.E = 1;
cParams.nu = 0;
cParams.bcVal = 1;
cParams.ndim = 2;
cParams.type = 'Elastic';
cParams.solverType='REDUCED';
cParams.solverMode = 'DISP';
cParams.solverCase = 'DIRECT';
cParams.scale = 'MACRO';

tolerance = 1e-10;
type = 'FORCE'; %"INTERNAL OR EXTERNAL"
%type = "INTERNAL";
TEST = TestingContinuumDamage(cParams,tolerance,type);
TEST.compareWithElasticProblem(cParams);