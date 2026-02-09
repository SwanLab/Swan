clc,clear,close all

%% INITIAL CONDITIONS FOR CONTINUING AN UNFINISHED ANALYSIS
% load("SENshearAT2_v2.mat")
% s.initialGuess.u = outputData.displacement.field;
% s.initialGuess.phi = outputData.damage.field.fun;

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

s.benchmark.mesh.type   = 'Rectangle';%'SENshear';
s.benchmark.mesh.length = 200;
s.benchmark.mesh.width  = 10;
s.benchmark.mesh.lN     = 200;
s.benchmark.mesh.wN     = 10;
s.benchmark.bc.u.type   = 'DisplacementTractionX';%'DisplacementShear';
s.benchmark.bc.u.values =  [0:1e-3:0.03];%[0:1e-5:0.02];
s.benchmark.bc.phi.type = 'DamageFixedLimitsX';%'DamageFixedLimitsX'; %DamageFree

s.matInfo.matType = 'Homogenized'; %'Analytic','Homogenized'
s.matInfo.degradationType = 'PhaseField'; %'PhaseField','SIMPALL'
s.matInfo.degradationSubType = 'AT'; %'AT','ATSplit','Rational','General'
s.matInfo.fileName = 'HoneycombBenchmark02'; 
s.matInfo.young   = 3*1e4;
s.matInfo.poisson = 0.2;
s.matInfo.Gc      = 0.008;
s.matInfo.sigmaMax =3;
s.l0 = 10; %10/5
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
save("UniaxialHoneyNew.mat",'outputData') %ACTIVATE TO SAVE DATA!
PhaseFieldPlotter(outputData);


s.matInfo.fileName = 'HexagonBenchmark02'; 
tester = TestingPhaseField(s);
outputData = tester.compute();
outputData.inputParameters = s;
save("UniaxialHexaNew.mat",'outputData')