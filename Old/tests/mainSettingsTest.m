
clc; clear variables; close all;
addpath(genpath(fileparts(mfilename('fullpath'))));

clear variables;

%% Main
filenames={
    'paramsTopOptProblem'
    };


for icases=1:size(filenames,1)
    clearvars -except filenames icases;
    close all;
    
    settings = SettingsTopOptProblem(filenames{icases});
    
    topOptProblem = TopOpt_Problem(settings);
    topOptProblem.computeVariables;
    topOptProblem.postProcess;
end
close all
