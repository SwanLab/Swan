clc,clear,close all

%% INITIAL CONDITIONS FOR CONTINUING AN UNFINISHED ANALYSIS
% load("SENshearAT2_v2.mat")
% s.initialGuess.u = outputData.displacement.field;
% s.initialGuess.phi = outputData.damage.field.fun;

%% GENERAL SETTINGS

s.monitoring.set =true;
s.monitoring.type = 'full'; %'reduced'
s.monitoring.print = true;

s.tolerance.u = 1e-12;
s.tolerance.phi = 1e-12;
s.tolerance.stag = 1e-12;
s.maxIter.u = 100;
s.maxIter.phi = 300;
s.maxIter.stag = 300;

s.benchmark.mesh.type   = 'Rectangle';%'SENshear';
s.benchmark.mesh.length= 1;
s.benchmark.mesh.width=0.1;
s.benchmark.mesh.lN= 50;
s.benchmark.mesh.wN= 5;
s.benchmark.bc.u.type   = 'DisplacementTractionXClamped';%'DisplacementShear';
s.benchmark.bc.u.values =  [0:1e-1:2];
s.benchmark.bc.phi.type = 'DamageFixedLimitsX';

s.matInfo.matType = 'Analytic'; %'Analytic','Homogenized'
s.matInfo.degradationType = 'PhaseField'; %'PhaseField','SIMPALL'
s.matInfo.degradationSubType = 'AT'; %'AT','AT2linear','General','IsoHomog'
s.matInfo.fileName = 'HexagonBenchmark03'; % (Only for Homogenized material)

s.matInfo.young   = 1;
s.matInfo.poisson = 0.3;
s.matInfo.Gc      = 1;
s.l0              = 0.1;


% s.matInfo.params.coeffs = [(4/pi)*(s.matInfo.Gc*s.matInfo.young)/(s.matInfo.sigmaMax^2 * s.l0), -0.5]; %(4/pi)
% s.matInfo.params.exp = 2;
% s.matInfo.params.fileName = 'HexagonBenchmark02'; % (Only for IsoHomog material)

s.dissipInfo.type = 'AT';
s.dissipInfo.constant = 2; % 2 AT2 / 8/3 AT1 / pi Wu 
s.dissipInfo.pExp = 2;
s.dissipInfo.xi = 0; % 0 AT2 / 1 AT1 / 2 Mix in type FullQuadratic
s.solver.type = 'Newton';
s.solver.tau  = 150;

sigmaPeak = sqrt((2*s.matInfo.Gc*s.matInfo.young)/(s.dissipInfo.constant *s.l0))

L = 1;
lch = (2*s.matInfo.Gc*s.matInfo.young)/(sigmaPeak^2)
lstr = L/lch
lreg = s.l0/lch

%% RUN
tester = TestingPhaseField(s);
outputData = tester.compute();
outputData.inputParameters = s;

%% SAVE + PLOT
%save("1ElemRational1_v2.mat",'outputData') %ACTIVATE TO SAVE DATA!
PhaseFieldPlotter(outputData);
% save('1ElemHexagonSimplified')
