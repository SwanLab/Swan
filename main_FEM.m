clc
clear all
addpath('.\FEM\');
%% Steps
% 1 - Run 'Main.m'
% 2 - Create object --> obj = Physical_Problem();
% 3 - Preprocess    --> obj.preProcess(filename_in);
% 4 - Compute       --> obj.computeVariables();
% 5 - Postprocess   --> obj.postProcess(filename_out);

%% New Commit:
% - New property attributes: more hermetic
% - New test: quadratic triangle

%% Main.m

% name_in  = 'CantileverToy_Tetrahedra';
% name_in  = 'CantileverToy_Hexahedra';
name_in = 'test2d_triangle';

name_out=strcat('results-',name_in);

folder_in  = 'Input';
folder_out = 'Output';

dir_in  = fullfile(pwd,folder_in);
dir_out = fullfile(pwd,folder_out);

filename_in  = fullfile(dir_in,name_in);
filename_out = fullfile(dir_out,name_out); 



obj = Physical_Problem();
% props.kappa = 1; props.mu = 0.4;
% obj.setMatProps(props);
obj.preProcess(filename_in);
obj.computeVariables();
obj.postProcess(filename_out);



