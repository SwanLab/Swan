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
s.monitoring.print = false;

s.tolerance.u = 1e-13;
s.tolerance.phi = 1e-6;
s.tolerance.stag = 1e-6;

s.benchmark.N = 10;
s.benchmark.type.mesh = '1Elem';
s.benchmark.type.bc = 'forceTraction';
s.benchmark.bcValues = [0:0.025:1.075];%[0:0.0001:0.1];

s.matInfo.matType = 'Analytic';
s.matInfo.degradationType = 'AT';
s.matInfo.fileName = 'CircleMicroDamagePerimeter'; 
s.matInfo.young   = 210;
s.matInfo.poisson = 0.3;
s.matInfo.Gc = 5e-3;
s.l0 = 0.1;

s.dissipInfo.type = 'PhaseFieldDissipationAT';
s.dissipInfo.pExp = 2;
s.solverType = 'Newton';

%% RUN
tester = TestingPhaseField(s);
outputData = tester.compute();
outputData.inputParameters = s;

%% SAVE + PLOT
%save("~/Documents/GitHub/Swan/SENshear_SquareArea.mat","outputData") %ACTIVATE TO SAVE DATA!
PhaseFieldPlotter(outputData);