clc,clear,close all

%% INITIAL CONDITIONS FOR CONTINUING AN UNFINISHED ANALYSIS
% load("nameOfOutputDataFile.mat")
% s.initialGuess.u = outputData.displacement.field;
% s.initialGuess.phi = outputData.damage.field;

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

s.benchmark.mesh.type   = 'SENtraction';
s.benchmark.mesh.length = 200;
s.benchmark.mesh.width  = 10;
s.benchmark.mesh.lN     = 200;
s.benchmark.mesh.wN     = 10;
s.benchmark.bc.u.type   = 'DisplacementTractionY';
s.benchmark.bc.u.values =  [0:1e-5:0.005,0.005:1e-6:0.1];
s.benchmark.bc.phi.type = 'DamageFree';%'DamageFixedLimitsX'; %DamageFree

s.matInfo.matType = 'Homogenized'; %'Analytic','Homogenized'
s.matInfo.degradationType = 'PhaseField'; %'PhaseField','SIMPALL'
s.matInfo.degradationSubType = 'General'; %'AT','ATSplit','Rational','General'
s.matInfo.fileName = 'HoneycombAreaNew'; 
s.matInfo.young   = 210; %3*1e4
s.matInfo.poisson = 0.3; %0.2
s.matInfo.Gc      = 2.7e-3; %0.008
s.matInfo.sigmaMax = 1.984; %3
s.l0 = 0.01; %10/5
s.matInfo.params.coeffs = [(4/pi)*(s.matInfo.Gc*s.matInfo.young)/(s.matInfo.sigmaMax^2 * s.l0), -0.5]; %(4/pi)
s.matInfo.params.exp = 2;

s.dissipInfo.type = 'AT';
s.dissipInfo.constant = 8/3; % 2 AT2 / 8/3 AT1 / pi Wu 
s.dissipInfo.pExp = 1;
s.dissipInfo.xi = 1; % 0 AT2 / 1 AT1 / 2 Mix in type FullQuadratic
s.solver.type = 'Gradient';
s.solver.tau  = 150;


%% RUN
tester = TestingPhaseField(s);
outputData = tester.compute();
outputData.inputParameters = s;

%% SAVE + PLOT
save("SENtractionHoney",'outputData') %ACTIVATE TO SAVE DATA!
PhaseFieldPlotter(outputData);

% Check reactions, and initial damage