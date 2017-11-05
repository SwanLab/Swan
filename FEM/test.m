clear all; close all; clc
%% TEST
% -

% Load the results for 2-d and 3-d tests
tests={'test2d_triangle';
    'test2d_quad';
    'test3d_hexahedra';
    'test3d_tetrahedra'};
% Parent directory
[parentdir,~,~] = fileparts(mfilename('fullpath'));

% Run Main.m
for i=1:length(tests)
    file_name = tests{i};
    file_name_in = strcat('./Input/',tests{i});
    
    load_file = strcat('./tests/',file_name);
    load(load_file)
    obj = Physical_Problem(file_name_in);
    obj.preProcess;
    obj.computeVariables;
    if sum(abs(obj.variables.d_u - d_u)) < 1e-6
        disp(strcat(file_name,' PASSED'));
    else
        disp(strcat(file_name,' FAILED'));
    end
end
