clc; clear variables; close all;
addpath(genpath(fileparts(mfilename('fullpath'))));

%% Test
run('test_fem.m');
run('test_topopt.m');
clear variables;

%% Settings
settings = Settings('CantileverTriangle_Case_4_1_3');

%% Main
tic
test = TopOpt_Problem(settings);
test.preProcess;
test.computeVariables;
toc
% test.postProcess;