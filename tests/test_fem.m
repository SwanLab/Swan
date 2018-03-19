clear all; close all; clc
%% TEST
%
% Material Lame Parameters
% kappa = .9107
% mu = .3446
%

% Load the results for 2-d and 3-d tests
tests_fem={'test2d_triangle';
    'test2d_quad';
    'test3d_hexahedra';
    'test3d_tetrahedra';
    'test2d_triangle_neo'};
% Parent directory
[parentdir,~,~] = fileparts(mfilename('fullpath'));


%% FEM TESTS
for i=1:length(tests_fem)
    file_name = tests_fem{i};
    file_name_in = strcat('./Input/',tests_fem{i});
    
    load_file = strcat('./tests/',file_name);
    load(load_file)
    obj = Physical_Problem(file_name);
    obj.preProcess;
    obj.computeVariables;
    if sum(abs(norm(obj.variables.d_u - d_u)/norm(d_u))) < 1e-6
        disp(strcat(file_name,' PASSED'));
    else
        disp(strcat(file_name,' FAILED'));
    end
end