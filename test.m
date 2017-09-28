clear all; close all; clc
%% TEST
% - 

% Load the results for 2-d and 3-dtest
tests={'test2d';
      'test3d'};
% Parent directory
[parentdir,~,~] = fileparts(pwd);

% Run Main.m
for i=1:length(tests)
file_name=tests{i};
load_file=strcat('./tests/',file_name);
load(load_file)
obj = MainFunc(file_name);

if sum(abs(obj.variables.displacement - d_u)) < 1e-6
    disp(strcat(file_name,' PASSED'));
else
    disp(strcat(file_name,' FAILED'));
end
end