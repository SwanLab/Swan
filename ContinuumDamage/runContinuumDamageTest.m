clc;clear;

load('ContinuumDamageTestInput.mat')

TEST = TestingContinuumDamage(s);
TEST.compareWithElasticProblem(s);