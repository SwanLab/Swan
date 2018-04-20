clc; clear variables; close all;
addpath(genpath(fileparts(mfilename('fullpath'))));

%% Test
% run('test_fem.m');
% run('test_topopt.m');
clear variables;
%% Settings

% settings = Settings('CantileverTriangle_Case_2_2_3');
settings=Settings('Vfrac07n5');%PruebaVfrac07P03   Vfrac05n5

%% Main
tic
test = TopOpt_Problem(settings);
test.preProcess;
test.computeVariables;
toc
test.postProcess;