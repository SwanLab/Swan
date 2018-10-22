%% INTEGRATION TEST =======================================================

clear; close all;

fprintf('Running TopOpt tests...\n')

%% Test Declaration -------------------------------------------------------

% tests_integration = {'test_circle_triangle','test_circle_quadrilateral','test_sphere_tetrahedra','test_sphere_hexahedra'};
tests_integration = {'test_cylinder_tetrahedra','test_cylinder_hexahedra'};

%% Run Integration Opt Tests ----------------------------------------------
for i = 1:length(tests_integration)
    clearvars -except tests_integration i
    %     tic
    file_name = tests_integration{i};
    file_name_in = strcat('./Input/',tests_integration{i});
    settings = Settings(file_name_in);
    load_file = strcat('./tests/',file_name);
    load(load_file)
    
    obj = TopOpt_Problem(settings);
    obj.preProcess;
    x = obj.x;
    
    filter_boundary = Filter_Boundary.create(obj.settings);
    filter_boundary.setupFromGiDFile(obj.settings.filename,obj.settings.ptype);
    filter_boundary.preProcess;
    A = filter_boundary.computeSurface(x);
    errorSurf = A/A0 - 1;
    if abs(errorSurf) < 4e-2
        cprintf('green',strcat(file_name,' PASSED.  Surface Error: ',num2str(errorSurf),'\n'));
    else
        cprintf('err',strcat(file_name,' FAILED. Surface Error: ',num2str(errorSurf),'\n'));
    end
    
    filter = Filter.create(obj.settings);
    filter.setupFromGiDFile(obj.settings.filename,obj.settings.ptype);
    filter.preProcess;
    V = filter.computeVolume(x);
    errorVol= V/V0 - 1;
    if abs(errorVol) < 6e-2
        cprintf('green',strcat(file_name,' PASSED.  Volume Error: ',num2str(errorVol),'\n'));
    else
        cprintf('err',strcat(file_name,' FAILED. Volume Error: ',num2str(errorVol),'\n'));
    end
    
    %     toc
    clear settings
end

fprintf('\nTopOpt tests completed.\n')
fprintf('\n-------------------------------------------\n\n')