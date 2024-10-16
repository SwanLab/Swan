clc;clear;

load('ContinuumDamageTestInput.mat')
tolerance = 1e-12;
type = "EXTERNAL"; %"INTERNAL OR EXTERNAL"
%type = "INTERNAL";
TEST = TestingContinuumDamage(s,tolerance,type);
TEST.compareWithElasticProblem(s);