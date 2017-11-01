clc
clear all
addpath('.\FEM\');
%% Steps
% 1 - Run 'Main.m'
% 2 - Create object --> obj = Physical_Problem();
% 3 - Preprocess    --> obj.preProcess(filename_in);
% 4 - Compute       --> obj.computeVariables();
% 5 - Postprocess   --> obj.postProcess(filename_out);

%% test
run('test.m')
clear variables
%% Main.m

triangle = Physical_Problem('CantileverToy_Triangular');
% props.kappa = 1; props.mu = 0.4;
% obj.setMatProps(props);
triangle.preProcess;
triangle.computeVariables;
triangle.postProcess;

tetrahedra=Physical_Problem('CantileverToy_Tetrahedra');
tetrahedra.preProcess;
tetrahedra.computeVariables;
tetrahedra.postProcess;


