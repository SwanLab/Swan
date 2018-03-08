clc; clear variables; close all;
addpath(genpath(fileparts(mfilename('fullpath'))));
%% Test
%run('test_fem.m');
run('test_topopt.m');
clear variables;
%% Settings
settings=Settings('CantileverTriangle_Case_1_2_1');
settings.filename='Bridge';
settings.target_parameters.Vfrac=0.2;
%% main
tic
test = TopOpt_Problem(settings);
test.preProcess;
test.computeVariables;
toc
test.postProcess;