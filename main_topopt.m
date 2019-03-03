
clc; clear variables; close all;
addpath(genpath(fileparts(mfilename('fullpath'))));

%% Test
% run('PlottingTests.m')
% run('FemTests.m');
% run('TopOptTests.m');
% run('UnfittedIntegrationTests.m')
% run('PrintingTests.m')
% run('AllTests.m')
clear variables;

%% Main
filenames={
    'CantileverTriangleCoarse_Case_1_1_1'
%     'ImprovedBridgeSYM_Case_5_1_4'     
%     'ThroneTetrahedraSYM_Case_1_1_5'
%     'ThroneTetrahedraSYM_Case_5_1_2'
%     'BikeTriangle_1_1_1'
    };


for icases=1:size(filenames,1)
    clearvars -except filenames icases;
    close all;
    settings = Settings(filenames{icases});
    
    topOptProblem = TopOpt_Problem(settings);
    topOptProblem.preProcess;
    topOptProblem.computeVariables;
    topOptProblem.postProcess;
end
close all
