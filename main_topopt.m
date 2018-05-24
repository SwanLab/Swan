clc; clear variables; close all;
addpath(genpath(fileparts(mfilename('fullpath'))));

%% Test
%run('test_fem.m');
%run('test_topopt.m');
clear variables;
%% Settings
%settings = Settings('CantileverTriangle_Case_2_2_1');
%settings = Settings('CantileverQuadrilateral_Case_1_1_1');
settings = Settings('CantileverTetrahedraFine_Case_1_2_1');
%settings = Settings('CantileverTetrahedra_Case_1_1_1');

%settings = Settings('Bulk1_3D');
%settings = Settings('BridgeTetrahedra_Case_2_1_1');
%settings=Settings('HoneyComb');%Vfrac07   Vfrac05n5 Bulk1 MicroTriangle_Case_3_7_1
%settings.printing=false;

%% Main
tic
test = TopOpt_Problem(settings);
test.preProcess;
test.computeVariables;
toc
test.postProcess;
