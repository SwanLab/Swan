clc,clear,close all

%% INITIAL CONDITIONS FOR CONTINUING AN UNFINISHED ANALYSIS
% load("nameOfOutputDataFile.mat")
% s.initialGuess.u = outputData.displacement.field;
% s.initialGuess.phi = outputData.damage.field;

%% LOAD PREDEFINED CASE
%load('case1ElemAT2.mat')

%% GENERAL SETTINGS
s.monitoring.set = true;
s.monitoring.type = 'full'; %'reduced'
s.monitoring.print = true;

s.tolerance.u = 1e-13;
s.tolerance.phi = 1e-6;
s.tolerance.stag = 1e-6;

s.benchmark.N = 10;
s.benchmark.type.mesh = 'SENshear';
s.benchmark.type.bc = 'displacementShear';
s.benchmark.bcValues = [1e-10:1e-5:5e-3, 5e-3:1e-6:1];
s.matInfo.E  = 210;
s.matInfo.nu = 0.3;
s.matInfo.Gc = 2.7e-3;
s.l0 = 0.01;

s.matInfo.matType = 'PhaseFieldHomog';
s.matInfo.fileName = 'CircleMicroDamagePerimeter'; 
s.matInfo.degradation = 'PhaseFieldDegradation';
s.dissipInfo.type = 'PhaseFieldDissipationAT';
s.dissipInfo.pExp = 2;
s.solverType = 'Gradient';

%% RUN
tester = TestingPhaseField(s);
outputData = tester.compute();
outputData.inputParameters = s;
save('caseSENshearCirclePerimeter','outputData')
%% SAVE + PLOT
% save("/home/gerard/Documents/GitHub/Swan/PhaseFieldExperiments/ResultsPaper1/SENshear/" + ...
%       "SENshear_SquareArea.mat","outputData") %ACTIVATE TO SAVE DATA!
PhaseFieldPlotter(outputData);