% function Main
clc; clear all;

%% File name
name_in  = 'toyExample.msh';
name_out = 'results_toyExample';

% name_in  = 'Cantileverbeam.msh';
% name_out = 'results_Cantileverbeam'; 

folder_in  = 'Input';
folder_out = 'Output';

dir_in  = fullfile(pwd,folder_in);
dir_out = fullfile(pwd,folder_out);

filename_in  = fullfile(dir_in,name_in);
filename_out = fullfile(dir_out,name_out); 

%% Physical Problem Object
% cantilever = Physical_Problem();
% cantilever.preProcess(filename_in);
toyExample = Physical_Problem();
toyExample.preProcess(filename_in);
toyExample.computeVariables();
toyExample.postProcess(filename_out);