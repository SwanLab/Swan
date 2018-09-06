clc; clear variables; close all;
addpath(genpath(fileparts(mfilename('fullpath'))));

%% Test
%run('test_fem.m');
run('test_topopt.m');
clear variables;
%% Settings
%settings = Settings('CantileverTriangle_Case_3_2_1');

%% Main
filenames ={ %Write several cases here to compute all
    'CantileverHexahedraCoarse_Case_1_1_1';
   %'CantileverQuadCoarse_Case_1_1_1';
 %   'ChairTetrahedraCoarse_Case_1_1_1';
       
  %'ImpBridgeHexahedra_Case_1_1_2';
 %'ImpCantileverHexahedra_Case_1_1_1';
 %'ImpCantileverHexahedra_Case_1_1_2';
 %'ImpCantileverHexahedra_Case_1_1_3';
};

for icases=1:size(filenames,1)
clearvars -except filenames icases iter;
close all;
settings=Settings(filenames{icases});

settings.nsteps=1;

test = TopOpt_Problem(settings);
test.preProcess;
test.computeVariables;
test.postProcess;
end