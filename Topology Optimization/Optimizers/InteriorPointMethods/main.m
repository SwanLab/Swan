clc; clear all; close all

% load bpopt libraries
addpath('bpopt');

problemNumber = 6;              % input number from 1 to 10
prob = ProblemCreator(problemNumber); 
prob.create();
bp = prob.x;
bp.mu = 100;

solveProblem = InteriorPointMethodsSolver(bp);
solveProblem.solve();
sol = solveProblem.sol;