clc; clear variables; close all;
addpath(genpath(fileparts(mfilename('fullpath'))));

%% Test
%run('test_fem.m');
%run('test_topopt.m');
clear variables;
%% Settings
%settings = Settings('CantileverTriangle_Case_3_2_1');

%% Main
filenames ={ %Write several cases here to compute all
   
%     'BridgeTriangleCoarse_Case_1_1_1';
};

for icases=1:size(filenames,1)
clearvars -except filenames icases;
close all;
settings=Settings(filenames{icases});
test = TopOpt_Problem(settings);
test.preProcess;
test.computeVariables;
test.postProcess;
end