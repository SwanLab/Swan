clc;clear;

load('inputTest.mat')
% load ('inputUnitMesh.mat');
% 
cParams.meshLength = 1;
cParams.meshWidth = 1;
cParams.meshN = 1;
cParams.meshM = 1;
cParams.E = 1;
cParams.nu = 0;
cParams.bcVal = 1;
cParams.ndim = 2;
cParams.type = 'Elastic';
cParams.solverType='REDUCED';
cParams.solverMode = 'DISP';
cParams.solverCase = 'DIRECT';
cParams.scale = 'MACRO';

tol = 1e-10;
% type = 'FORCE'; %"DISP OR FORCE"
type = 'FORCE';
TEST = TestingContinuumDamage(cParams,tol,type);
TEST.compareWithElasticProblem();