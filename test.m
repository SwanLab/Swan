clear; close all; clc
%% TEST
% - 
tic
% Load the results for 2-d and 3-d tests
tests={'test2d_triangle';
        'test2d_quad';
        'test3d_hexahedra';
      'test3d_tetrahedra';
      'test3d_tetrahedra_bis'};
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
toc