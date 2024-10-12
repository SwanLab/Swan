clc;clear;

load('ContinuumDamageTestInput.mat')
tolerance = 1e-12;
TEST = TestingContinuumDamage(s,tolerance);
TEST.compareWithElasticProblem(s);