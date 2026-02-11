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

s.benchmark.mesh.type   = 'SENtraction';%'SENshear';
s.benchmark.mesh.length = 200;
s.benchmark.mesh.width  = 10;
s.benchmark.mesh.lN     = 200;
s.benchmark.mesh.wN     = 10;
s.benchmark.bc.u.type   = 'DisplacementTractionYClamped';%'DisplacementShear';
s.benchmark.bc.u.values =  [0:1e-5:0.005,0.005:1e-6:0.0065];
s.benchmark.bc.phi.type = 'DamageFree';%'DamageFixedLimitsX'; %DamageFree

s.matInfo.matType = 'Analytic'; %'Analytic','Homogenized'
s.matInfo.degradationType = 'PhaseField'; %'PhaseField','SIMPALL'
s.matInfo.degradationSubType = 'AT'; %'AT','ATSplit','Rational','General'
s.matInfo.fileName = 'HoneycombBenchmark02'; 
s.matInfo.young   = 210;
s.matInfo.poisson = 0.3;
s.matInfo.Gc      = 2.7e-3;
s.matInfo.sigmaMax =2.44542;
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
save("SENtraction001AT1.mat",'outputData') %ACTIVATE TO SAVE DATA!
PhaseFieldPlotter(outputData);

s.matInfo.degradationSubType = 'AT2linear';
s.dissipInfo.constant = 2;
tester = TestingPhaseField(s);
outputData = tester.compute();
outputData.inputParameters = s;
save("SENtraction001AT2.mat",'outputData') 

s.matInfo.degradationSubType = 'General';
s.matInfo.sigmaMax =2.44542;
s.dissipInfo.constant = pi;
tester = TestingPhaseField(s);
outputData = tester.compute();
outputData.inputParameters = s;
save("SENtraction001Rational24.mat",'outputData') 

s.matInfo.sigmaMax = 5;
s.dissipInfo.constant = pi;
tester = TestingPhaseField(s);
outputData = tester.compute();
outputData.inputParameters = s;
save("SENtraction001Rational5.mat",'outputData') 