clc;clear;

load('ContinuumDamageTestInput.mat')
% 
% s.meshLength = 1;
% s.meshWidth = 1;
% s.meshN = 2;
% s.meshM = 2;
% s.E = 1;
% s.nu = 0;
% s.bcVal = 1;
% s.ndim = 2;
% s.type = 'Elastic';
% s.solverType='REDUCED';
% s.solverMode = 'DISP';
% s.solverCase = 'DIRECT';
% s.scale = 'MACRO';

tolerance = 1e-10;
type = 'FORCE'; %"DISP OR FORCE"

TEST = TestingContinuumDamage(s,tolerance,type);
TEST.compareWithElasticProblem();