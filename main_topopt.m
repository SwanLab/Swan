clc; clear variables; close all;
addpath(genpath(fileparts(mfilename('fullpath'))));

%% Test
% run('test_fem.m');
% run('test_topopt.m');
clear variables;

%% Settings
settings = Settings('CantileverTriangle_Case_3_2_2');
settings.printing = 1;
settings.plotting = 0;
settings.monitoring = 0;

%% Main
tic
test = TopOpt_Problem(settings);
test.preProcess;
test.computeVariables;
toc
test.postProcess;