clc; clear variables; close all;
addpath(genpath(fileparts(mfilename('fullpath'))));

%% Test
run('test_fem.m');
run('test_topopt.m');
clear variables;
%% Settings
settings = Settings('CantileverTriangle_Case_2_1_3');
%settings = Settings('CantileverQuadrilateral_Case_1_2_1');
%settings = Settings('CantileverTetrahedra_Case_2_1_2');
%% Main
tic
test = TopOpt_Problem(settings);
test.preProcess;
test.computeVariables;
toc
test.postProcess;