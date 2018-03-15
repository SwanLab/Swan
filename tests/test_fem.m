clear all; close all; clc
%% TEST
% -

% Load the results for 2-d and 3-d tests
tests_fem={'test2d_triangle';
    'test2d_quad';
    'test3d_hexahedra';
    'test3d_tetrahedra'
    'test2d_micro'};

%% FEM TESTS
for i=1:length(tests_fem)
    file_name = tests_fem{i};
    file_name_in = strcat('./Input/',tests_fem{i});
    
    load_file = strcat('./tests/',file_name);
    load(load_file);   
    
    obj = Physical_Problem(file_name);
    if obj.mesh.scale == 'MACRO'
        obj.preProcess;
        obj.computeVariables;
        if sum(abs(obj.variables.d_u - d_u)) < 1e-6
            cprintf('green',strcat(file_name,' PASSED\n'));
        else
             cprintf('err',strcat(file_name,' FAILED\n'));
        end
    else
        obj = Physical_Problem_Micro(file_name);
        obj.preProcess;
        obj.computeChomog;
        if sum(abs(obj.variables.Chomog- Chomog)) < 1e-6
             cprintf('green',strcat(file_name,' PASSED\n'));
        else
             cprintf('err',strcat(file_name,' FAILED\n'));
        end
    end
    
end