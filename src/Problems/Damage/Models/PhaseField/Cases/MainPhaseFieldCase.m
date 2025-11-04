clc,clear,close all

%% INITIAL CONDITIONS FOR CONTINUING AN UNFINISHED ANALYSIS
% load("nameOfOutputDataFile.mat")
% s.initialGuess.u = outputData.displacement.field;
% s.initialGuess.phi = outputData.damage.field;

%% LOAD PREDEFINED CASE
%load('caseSENtractionAT2.mat')

%% GENERAL SETTINGS
s = []; 

s.monitoring.set = true;
s.monitoring.type = 'full'; %'reduced'
s.monitoring.print = true;

s.tolerance.u = 1e-6;
s.tolerance.phi = 1e-10;
s.tolerance.stag = 1e-8;
s.maxIter.u = 100;
s.maxIter.phi = 100;
s.maxIter.stag = 300;

s.benchmark.mesh.type = '1Elem';
s.benchmark.bc.type   = 'DisplacementTractionX';
s.benchmark.bc.values =  [0:1e-4:0.1];

s.matInfo.matType = 'Analytic'; %'Analytic','Homogenized'
s.matInfo.degradationType = 'PhaseField'; %'PhaseField','SIMPALL'
s.matInfo.degradationSubType = 'Rational'; %'AT','ATSplit',,'Rational','General'
s.matInfo.fileName = 'L0_variation/DegSqr15L0LHS'; 
s.matInfo.young   = 210;
s.matInfo.poisson = 0.3;
s.matInfo.Gc      = 5e-3;
s.matInfo.sigmaMax = 1.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E = s.matInfo.young;
nu = s.matInfo.poisson;
Gc = s.matInfo.Gc;
sigma = s.matInfo.sigmaMax;

lch   = (2*E*Gc)/(sigma^2);

k  = 210./(2.*(1-0.3));
mu = 210./(2.*(1+0.3));
slope = - k - k^2/mu - mu - (2*mu^2 + k)/k;
lhs = -2*(3/8)*(5e-3/slope)*(210/sigma^2);


% s.matInfo.Gc = 5e-3*(sigma^2/sigmach^2);
% s.l0 = 0.1;

s.l0 = lch;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


s.dissipInfo.type = 'PhaseFieldDissipationAT';
s.dissipInfo.pExp = 1;
s.solver.type = 'Gradient';
s.solver.tau  = 150;


%% RUN
tester = TestingPhaseField(s);
outputData = tester.compute();
outputData.inputParameters = s;

%% SAVE + PLOT
save("Sqr15L02Lch",'outputData') %ACTIVATE TO SAVE DATA!
PhaseFieldPlotter(outputData);
