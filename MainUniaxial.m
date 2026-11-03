clc,clear,close all

%% INITIAL CONDITIONS FOR CONTINUING AN UNFINISHED ANALYSIS
% load("SENshearAT2_v2.mat")
% s.initialGuess.u = outputData.displacement.field;
% s.initialGuess.phi = outputData.damage.field.fun;

%% GENERAL SETTINGS

s.monitoring.set =true;
s.monitoring.type = 'full'; %'reduced'
s.monitoring.print = true;

s.tolerance.u = 1e-6;
s.tolerance.phi = 1e-6;
s.tolerance.stag = 1e-6;
s.maxIter.u = 100;
s.maxIter.phi = 300;
s.maxIter.stag = 300;

s.benchmark.mesh.type   = 'Rectangle';%'SENshear';
s.benchmark.mesh.length= 1;
s.benchmark.mesh.width=0.1;
s.benchmark.mesh.lN= 50;
s.benchmark.mesh.wN= 5;

s.benchmark.bc.u.type   = 'DisplacementTractionXClamped';%'DisplacementShear';
s.benchmark.bc.u.values =  [0:0.1:0.72, 0.72:1e-2:1];
s.benchmark.bc.phi.type = 'DamageFixedLimitsX';

s.matInfo.matType = 'Analytic'; %'Analytic','Homogenized'
s.matInfo.degradationType = 'PhaseField'; %'PhaseField','SIMPALL'
s.matInfo.degradationSubType = 'Rational'; %'AT','AT2linear','General','IsoHomog'

sigmaPeak = 1;
s.matInfo.young   = 2;
s.matInfo.poisson = 0.3;
s.matInfo.Gc      = 1;
s.l0              = 0.5;
s.matInfo.params.gamma = (2*s.matInfo.Gc*s.matInfo.young)/(pi*s.l0*sigmaPeak^2)
s.matInfo.params.gamma = 1.2
sPeak = sqrt((2*s.matInfo.Gc*s.matInfo.young)/(pi*s.l0*s.matInfo.params.gamma))


L = 1;
lch = pi*s.matInfo.params.gamma*s.l0
lch = (2*s.matInfo.Gc*s.matInfo.young)/(sigmaPeak^2)
lstr = L/lch
lreg = s.l0/lch

s.dissipInfo.type = 'FullQuadratic';
s.dissipInfo.xi = 2;
s.solver.type = 'Gradient';
s.solver.tau  = 150;


%% RUN
tester = TestingPhaseField(s);
outputData = tester.compute();
outputData.inputParameters = s;

%% SAVE + PLOT
%save("1ElemRational1_v2.mat",'outputData') %ACTIVATE TO SAVE DATA!
PhaseFieldPlotter(outputData);
