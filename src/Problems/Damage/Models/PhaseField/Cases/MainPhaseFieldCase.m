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
s.maxIter.phi = 300;
s.maxIter.stag = 300;

s.benchmark.mesh.type   = '1Elem';
% s.benchmark.mesh.length = 200;
% s.benchmark.mesh.width  = 10;
% s.benchmark.mesh.lN     = 200;
% s.benchmark.mesh.wN     = 10;

s.benchmark.bc.phi.type = 'DamageFree';
s.benchmark.bc.phi.isAdaptive = false;

s.benchmark.bc.u.type   = 'DisplacementTractionY';
s.benchmark.bc.u.isAdaptive = false;
s.benchmark.bc.u.values =  [0:5e-3:0.1];
% s.benchmark.bc.u.isAdaptive = true;
% s.benchmark.bc.u.initialValue = 0;
% s.benchmark.bc.u.finalValue = 0.03;

s.matInfo.matType = 'Homogenized'; %'Analytic','Homogenized'
s.matInfo.degradationType = 'PhaseField'; %'PhaseField','SIMPALL'
s.matInfo.degradationSubType = 'AT'; %'AT','AT2linear','General'
s.matInfo.fileName = 'HorizontalCrack'; 
s.matInfo.young   = 210;
s.matInfo.poisson = 0.3;
s.matInfo.Gc      = 5e-3;
s.matInfo.sigmaMax = 1.984;
s.l0 = 0.1;
s.matInfo.params.coeffs = [(4/pi)*(s.matInfo.Gc*s.matInfo.young)/(s.matInfo.sigmaMax^2 * s.l0), -0.5]; %(4/pi)
s.matInfo.params.exp = 2;

s.dissipInfo.type = 'AT';
s.dissipInfo.constant = 8/3; % 2 AT2 / 8/3 AT1 / pi Rational 
s.dissipInfo.pExp = 1;
s.dissipInfo.xi = 2; % 0 AT2 / 1 AT1 / 2 Mix in type FullQuadratic
s.solver.type = 'Gradient';
s.solver.tau  = 150;


%% Orientation plots
angleInfo.alpha = [];
angleInfo.beta  = [];
angleInfo.theta = [];
angleInfo.chi   = [];
angleInfo.phi   = [];
angleInfo.phiOld= [];
save('angleInfo','angleInfo')

%% RUN
tester = TestingPhaseField(s);
outputData = tester.compute();
outputData.inputParameters = s;

%% SAVE + PLOT
%save("1ElemRational1_v2.mat",'outputData') %ACTIVATE TO SAVE DATA!
PhaseFieldPlotter(outputData);

angleInfo = outputData.angleInfo;

for i=1:4
figure(100+i)
hold on
yyaxis right
p1 = plot(angleInfo.phi(i,:),'Color','Y');
p2 = plot(angleInfo.phiOld(i,:));
ylim([-0.1 1.1])
legend('phi')
yyaxis left
p3 = plot(rad2deg(angleInfo.alpha(i,:)),'Color','R');
p4 = plot(rad2deg(angleInfo.beta(i,:)),'Color','B');
p5 = plot(rad2deg(angleInfo.theta(i,:)),'Color','G');
p6 = plot(rad2deg(angleInfo.chi(i,:)),'Color','K');
ylim([-5 365])
legend([p1,p2,p3,p4,p5,p6],'phi','phiOld','alpha','beta','theta','chi')
end

