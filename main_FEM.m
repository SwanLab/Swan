clc
clear variables
addpath(genpath('.\FEM\'));
addpath(genpath('.\Input\'));
%% Steps
% 1 - Run 'Main.m'
% 2 - Create object --> obj = Physical_Problem(filename);
% 3 - Preprocess    --> obj.preProcess;
% 4 - Compute       --> obj.computeVariables;
% 5 - Postprocess   --> obj.postProcess;
%% test
run('test.m')
clear variables
%% Main.m
name_in = 'CantileverToy_Triangular';
name_out = strcat('results-',name_in);

triangle = Physical_Problem;
% props.kappa = 1; props.mu = 0.4;
% obj.setMatProps(props);
triangle.preProcess(name_in);
triangle.computeVariables;
triangle.postProcess(name_out);

name_in = 'CantileverToy_Tetrahedra';
name_out = strcat('results-',name_in);

tetrahedra = Physical_Problem;
tetrahedra.preProcess(name_in);
tetrahedra.computeVariables;
tetrahedra.postProcess(name_out);


