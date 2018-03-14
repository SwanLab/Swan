clc; clear variables; close all;
addpath(genpath(fileparts(mfilename('fullpath'))));
%% Test
run('test_fem.m');
run('test_topopt.m');
clear variables;
%% Settings
settings=Settings('CantileverTriangle_Case_1_1_3');
% settings=Settings('MicroTriangle_Case_1_6_1');
%% main
tic
test = TopOpt_Problem(settings);
test.preProcess;
test.computeVariables;
toc
% test.postProcess;