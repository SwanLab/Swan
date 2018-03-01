clc; clear variables; close all;
addpath(genpath(fileparts(mfilename('fullpath'))));
%% Test
%run('test_fem.m');
%run('test_topopt.m');
clear variables;
%% Settings
settings=Settings('Case4','CantileverBeam_Triangle_Linear_Fine');
%% main

tic
test = TopOpt_Problem(settings);
test.preProcess;
test.computeVariables;
toc
test.postProcess;


