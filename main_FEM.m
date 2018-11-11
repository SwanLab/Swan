clc
clear variables
% 
addpath(genpath(fileparts(mfilename('fullpath'))))
% 
% %% Steps
% % 1 - Run 'Main.m'
% % 2 - Create object --> obj = Elastic_Problem(filename);
% % 3 - Preprocess    --> obj.preProcess;
% % 4 - Compute       --> obj.computeVariables;
% % 5 - Postprocess   --> obj.postProcess;
% %% test
run('test_fem.m')
% % test
clear variables
%% Main.m

physProblem = Elastic_Problem('Cantilever_quad_coarse');

tic
physProblem.preProcess;
physProblem.computeVariables;
physProblem.print;
toc


