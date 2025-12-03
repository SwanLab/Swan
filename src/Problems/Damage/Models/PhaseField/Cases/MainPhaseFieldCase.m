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
s.tolerance.phi = 1e-6;
s.tolerance.stag = 1e-6;
s.maxIter.u = 100;
s.maxIter.phi = 100;
s.maxIter.stag = 300;

s.benchmark.mesh.type   = 'Xstructure';
s.benchmark.mesh.length = 1;
s.benchmark.mesh.width  = 1;
s.benchmark.mesh.lN     = 200;
s.benchmark.mesh.wN     = 200;
s.benchmark.bc.u.type   = 'DisplacementShear';
s.benchmark.bc.u.values =  [0:1e-2:1];
s.benchmark.bc.phi.type   = 'DamageSEN'; %DamageFree

s.matInfo.matType = 'Analytic'; %'Analytic','Homogenized'
s.matInfo.degradationType = 'PhaseField'; %'PhaseField','SIMPALL'
s.matInfo.degradationSubType = 'ATSplit'; %'AT','ATSplit',,'Rational','General'
s.matInfo.fileName = 'L0_variation/DegSqr15lHS'; 
s.matInfo.young   = 210;
s.matInfo.poisson = 0.3;
s.matInfo.Gc      = 5e-3;
s.matInfo.sigmaMax = 1.5;
s.l0 = 0.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


E = s.matInfo.young;
nu = s.matInfo.poisson;
Gc = s.matInfo.Gc;
sigma = s.matInfo.sigmaMax;

lch   = (2*E*Gc)/(sigma^2);

C11 = E/((1+nu)*(1-nu));
k  = 210./(2.*(1-0.3));
mu = 210./(2.*(1+0.3));
slope = (- k - k^2/mu - mu - (2*mu^2 + k)/k)/C11;
lhs = -2*(3/8)*(Gc/slope)*(E/sigma^2);


% s.matInfo.Gc = 5e-3*(sigma^2/sigmach^2);
% s.l0 = 0.1;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


s.dissipInfo.type = 'PhaseFieldDissipationAT';
s.dissipInfo.pExp = 2;
s.solver.type = 'Newton';
s.solver.tau  = 150;


%% RUN
tester = TestingPhaseField(s);
outputData = tester.compute();
outputData.inputParameters = s;

%% SAVE + PLOT
save("Xstructure",'outputData') %ACTIVATE TO SAVE DATA!
PhaseFieldPlotter(outputData);
print(outputData.damage.field.fun,'Damage')