clc;clear;

load('ContinuumDamageTestInput.mat')
tolerance = 1e-3;

TEST = TestingContinuumDamage(s,tolerance);
TEST.compareWithElasticProblem(s);