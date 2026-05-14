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
s.maxIter.u = 10;
s.maxIter.phi = 10;
s.maxIter.stag = 5;

s.benchmark.mesh.type   = 'Rectangle';
s.benchmark.mesh.length = 200;
s.benchmark.mesh.width  = 10;
s.benchmark.mesh.lN     = 200;
s.benchmark.mesh.wN     = 10;
s.benchmark.bc.phi.type = 'DamageFixedLimitsX';
s.benchmark.bc.phi.isAdaptive = false;
s.benchmark.bc.u.type   = 'DisplacementTractionXClamped';
% s.benchmark.bc.u.adaptive = false;
% s.benchmark.bc.u.values =  [0:1e-3:0.005,0.005:1e-6:0.0065];
s.benchmark.bc.u.isAdaptive = true;
s.benchmark.bc.u.initialValue = 0;
s.benchmark.bc.u.finalValue = 0.03;

s.matInfo.matType = 'Analytic'; %'Analytic','Homogenized'
s.matInfo.degradationType = 'PhaseField'; %'PhaseField','SIMPALL'
s.matInfo.degradationSubType = 'AT'; %'AT','AT2linear','General'
s.matInfo.fileName = 'HoneycombBenchmark02'; 
s.matInfo.young   = 3e4;
s.matInfo.poisson = 0.2;
s.matInfo.Gc      = 8e-3;
s.matInfo.sigmaMax = 3;
s.l0 = 10;
s.matInfo.params.coeffs = [(4/pi)*(s.matInfo.Gc*s.matInfo.young)/(s.matInfo.sigmaMax^2 * s.l0), -0.5]; %(4/pi)
s.matInfo.params.exp = 2;

s.dissipInfo.type = 'AT';
s.dissipInfo.constant = 8/3; % 2 AT2 / 8/3 AT1 / pi Rational 
s.dissipInfo.pExp = 1;
s.dissipInfo.xi = 2; % 0 AT2 / 1 AT1 / 2 Mix in type FullQuadratic
s.solver.type = 'Gradient';
s.solver.tau  = 150;

%% RUN
tester = TestingPhaseField(s);
outputData = tester.compute();
outputData.inputParameters = s;

%% SAVE + PLOT
%save("1ElemRational1_v2.mat",'outputData') %ACTIVATE TO SAVE DATA!
PhaseFieldPlotter(outputData);
