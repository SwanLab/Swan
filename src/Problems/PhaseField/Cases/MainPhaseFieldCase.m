clc,clear,close all

%% INITIAL CONDITIONS FOR CONTINUING AN UNFINISHED ANALYSIS
% load("nameOfOutputDataFile.mat")
% s.initialGuess.u = outputData.displacement.field;
% s.initialGuess.phi = outputData.damage.field;

%% LOAD PREDEFINED CASE
%load('caseSENtractionAT2.mat')

%% GENERAL SETTINGS
s.monitoring.set = true;
s.monitoring.type = 'full'; %'reduced'
s.monitoring.print = true;

s.tolerance.u = 1e-6;
s.tolerance.phi = 1e-6;
s.tolerance.stag = 1e-6;
s.maxIter.u = 100;
s.maxIter.phi = 100;
s.maxIter.stag = 300;

s.benchmark.N = 10;
s.benchmark.type.mesh = '1Elem';
s.benchmark.type.bc = 'displacementTraction';
s.benchmark.bcValues = [0:1e-3:1];

s.matInfo.matType = 'Analytic';
s.matInfo.degradationType = 'AT';
s.matInfo.fileName = 'SquarePerimeter'; 
s.matInfo.young   = 210;
s.matInfo.poisson = 0.3;
s.matInfo.Gc = 5e-3;
s.l0 = 0.1;

s.dissipInfo.type = 'PhaseFieldDissipationAT';
s.dissipInfo.pExp = 1;
s.solver.type = 'Newton';
s.solver.tau  = 150;


%% RUN
tester = TestingPhaseField(s);
outputData = tester.compute();
outputData.inputParameters = s;

%% SAVE + PLOT
% save("AT2shear_X",'outputData') %ACTIVATE TO SAVE DATA!
PhaseFieldPlotter(outputData);