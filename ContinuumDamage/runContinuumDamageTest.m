clc;clear;

s.meshLength = 10;
s.meshWidth = 1;
s.meshN = 10;
s.meshM = 10;
s.bcVal = 1;
s.ndim =2;
s.E = 10;
s.nu = 1/4;
s.type      = 'Elastic';
s.solverType = 'REDUCED';
s.solverMode = 'DISP';
s.solverCase = 'DIRECT';
s.scale = 'MACRO';

TEST = TestingContinuumDamage(s);
TEST.compareWithElasticProblem(s);