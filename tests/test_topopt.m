%% TOP OPT TEST ===========================================================

clear; close all;

% Test Declaration --------------------------------------------------------
tests_topopt={'test_cantilever','test_gripping'};

% Run Top Opt Tests -------------------------------------------------------
for i=1:length(tests_topopt)
    tic
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
        warning(strcat(file_name,' FAILED'));
    end
    toc
end