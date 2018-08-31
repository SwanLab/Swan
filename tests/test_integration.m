%% INTEGRATION TEST =======================================================

clear; close all;

fprintf('Running TopOpt tests...\n')

%% Test Declaration -------------------------------------------------------

tests_integration = {'test_sphere_tetrahedra','test_sphere_hexahedra'};

%% Run Integration Opt Tests ----------------------------------------------
for i = 1:length(tests_integration)
    clearvars -except tests_integration i
    tic
    file_name = tests_integration{i};
    file_name_in = strcat('./Input/',tests_integration{i});
    settings = Settings(file_name_in);
    load_file = strcat('./tests/',file_name);
    load(load_file)
    
    obj = TopOpt_Problem(settings);
    obj.preProcess;
    x = obj.x;
    
    filter =  Filter.create(obj.settings);
    filter.preProcess;
    A = filter.computeFacetSurface(x);
    A0 = 4*pi;
    A_A0 = A/A0;
    
    error = (A_A0 - A_A0_ref)/A_A0_ref;
    if error < 1e-9
        cprintf('green',strcat(file_name,' PASSED.  Error: ',num2str(error),'\n'));
    else
        cprintf('err',strcat(file_name,' FAILED. Error: ',num2str(error),'\n'));
    end
    toc
    clear settings
end

fprintf('\nTopOpt tests completed.\n')
fprintf('\n-------------------------------------------\n\n')