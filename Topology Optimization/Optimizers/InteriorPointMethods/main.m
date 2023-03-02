clc; clear all; close all

% load bpopt libraries
addpath('bpopt')

% Solve problem 1
prob = bp_create(10);  % create problem
prob.mu = 10;          % change the initial barrier term
sol = bpopt(prob);     % solve problem

% Solve problem 2
%prob = bp_create(6);   % create problem
%prob.mu = 100;         % change the initial barrier term
%sol = bpopt(prob);     % solve problem
