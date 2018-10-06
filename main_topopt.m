clc; clear variables; close all;
addpath(genpath(fileparts(mfilename('fullpath'))));

%% Test
run('test_fem.m');
run('test_topopt.m');
run('test_integration.m')
clear variables;

%% Main
filenames={
'CantileverbeamHexahedraSYM_Case_5_1_6'
    
% 'CantileverTetrahedraSYM_Case_1_1_1'
% 'CantileverTetrahedraSYM_Case_2_1_1'
% 'CantileverTetrahedraSYM_Case_3_1_1'
% 'CantileverTetrahedraSYM_Case_4_1_1'
% 'CantileverTetrahedraSYM_Case_5_1_1'

% 'CantileverTetrahedraSYM_Case_1_2_1'
% 'CantileverTetrahedraSYM_Case_2_2_1'
% 'CantileverTetrahedraSYM_Case_3_2_1'
% 'CantileverTetrahedraSYM_Case_4_2_1'
% 'CantileverTetrahedraSYM_Case_5_2_1'

% 'ImprovedBridgeSYM_Case_5_1_4'

%     'ThroneTetrahedraSYM_Case_1_1_5'
%     'ThroneTetrahedraSYM_Case_5_1_2'
};

for icases=1:size(filenames,1)
    tic
    clearvars -except filenames icases;
    close all;
    settings=Settings(filenames{icases});
    
%     try
        test = TopOpt_Problem(settings);
        test.preProcess;
        test.computeVariables;
        test.postProcess;
%     end
end
close all