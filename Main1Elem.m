clc,clear,close all

% GENERAL SETTINGS

s.monitoring.set = true;
s.monitoring.type = 'full'; %'reduced'
s.monitoring.print = true;

s.tolerance.u = 1e-6;
s.tolerance.phi = 1e-6;
s.tolerance.stag = 1e-6;
s.maxIter.u = 100;
s.maxIter.phi = 300;
s.maxIter.stag = 300;

s.benchmark.mesh.type   = '1Elem';
s.benchmark.bc.u.type   = 'DisplacementTractionX';
s.benchmark.bc.u.values =  [0:1e-4:0.1];
s.benchmark.bc.phi.type = 'DamageFree';

s.matInfo.matType = 'Analytic'; %'Analytic','Homogenized'
s.matInfo.degradationType = 'PhaseField'; %'PhaseField','SIMPALL'
s.matInfo.degradationSubType = 'General'; %'AT','ATSplit','Rational','General'
s.matInfo.fileName = 'HoneycombBenchmark02'; 
s.matInfo.young   = 210;
s.matInfo.poisson = 0.3;
s.matInfo.Gc      = 5e-3;
s.matInfo.sigmaMax = 1.984;
s.l0 = 0.1; %10/5
s.matInfo.params.coeffs = [(4/pi)*(s.matInfo.Gc*s.matInfo.young)/(s.matInfo.sigmaMax^2 * s.l0), -0.5]; %(4/pi)
s.matInfo.params.exp = 2;

s.dissipInfo.type = 'AT';
s.dissipInfo.constant = pi; % 2 AT2 / 8/3 AT1 / pi Rational 
s.dissipInfo.pExp = 1;
s.dissipInfo.xi = 2; % 0 AT2 / 1 AT1 / 2 Mix in type FullQuadratic
s.solver.type = 'Gradient';
s.solver.tau  = 150;


% RUN
tester = TestingPhaseField(s);
outputData = tester.compute();
outputData.inputParameters = s;

% SAVE + PLOT
save("1ElemRational_noUB.mat",'outputData') %ACTIVATE TO SAVE DATA!
PhaseFieldPlotter(outputData); 