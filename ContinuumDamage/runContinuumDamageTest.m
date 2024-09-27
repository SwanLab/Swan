clc;clear;

s.meshLength = 1;
s.meshWidth = 1;
s.meshN = 10;
s.meshM = 10;
s.bcVal = 1;
s.E = 10;
s.nu = 10;
s.type      = 'Elastic';
s.solverType = 'REDUCED'; %MIRAR COM INFLUEIX!!!
s.solverMode = 'DISP';
s.solverCase = 'DIRECT';
s.scale = 'MACRO';

TEST = TestingContinuumDamage(s);
TEST.compareWithElasticProblem(s);