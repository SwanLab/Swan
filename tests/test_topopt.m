clear all; close all;
%% TEST
tests_topopt={'test_gripping'};

%% TOP OPT TESTS
for i=1:length(tests_topopt)
    file_name = tests_topopt{i};
    file_name_in = strcat('./Input/',tests_topopt{i});
    load_file = strcat('./tests/',file_name);
    run(tests_topopt{i})
    load(load_file)
    settings.filename=file_name_in;
    obj = TopOpt_Problem(settings);
    obj.preProcess;
    obj.computeVariables;
    if sum(abs(obj.x - x)) < 1e-6
        disp(strcat(file_name,' PASSED'));
    else
        disp(strcat(file_name,' FAILED'));
    end
end