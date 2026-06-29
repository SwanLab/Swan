clear
clc
close all

%% Material

inputData.young   = 120e9;
inputData.poisson = 0.3;

%% Cohesive law

inputData.lawType = 'TractionBiliniarCoupled';

% Parameters that remain constant
inputData.tau0Shear        = 30e6;
inputData.firstCritEnergy  = 260;
inputData.secondCritEnergy = 1002;
inputData.eta              = 2;

% Alternative parameters
inputData.jumpCrit  = 1.25e-7;
inputData.jumpFinal = 0.025e-3;

inputData.Kcoh = 1e13;

%% Boundary conditions

inputData.bcValues = [ ...
    0      :1e-4   :1.15e-3,...
    1.15e-3   :2e-6  :1.5e-3,...
    1.5e-3:1e-5   :2.5e-3];

%% Problem definition

inputData.problemType = 'DoubleCantileverBeam';
inputData.solverType  = 'Newton';

%% Geometry

inputData.l = 150e-3;
inputData.h = 1.55*2e-3;

inputData.yCohLine    = 0.5*inputData.h;
inputData.xCohLineMax = 115e-3;

%% Mesh
    inputData.nx = 1500;
    inputData.ny = 10;

%% Cases to analyze

% tau0NormalValues = [30e6, 60e6, 15e6];

tau0NormalValues = [60e6];

results = struct();

%% Run simulations

for iCase = 1:length(tau0NormalValues)

    inputData.tau0Normal = tau0NormalValues(iCase);

    fprintf('\n=====================================\n');
    fprintf('Running case %d/%d: tau0Normal = %.1f MPa\n', ...
        iCase, length(tau0NormalValues), inputData.tau0Normal/1e6);
    fprintf('=====================================\n');

    problem = IterativeTutorialCohesive(inputData);

    results(iCase).tau0Normal = inputData.tau0Normal;
    results(iCase).problem    = problem;
    results(iCase).output     = problem.output;

end