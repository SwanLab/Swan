%% RunLastCaseAndPlotComparison
% Carrega os tres casos concluidos, executa somente o TwoVariable
% e gera novamente os graficos e a tabela comparativa.

clear;
close all;
clc;

clear classes;
rehash;

outputFolder = fullfile(pwd,'TopologyComparisonResults');

if ~exist(outputFolder,'dir')
    error('Pasta nao encontrada: %s',outputFolder);
end

fileSIMP = fullfile(outputFolder,'Result_SIMP.mat');
fileB    = fullfile(outputFolder,'Result_BOnly.mat');
fileRho  = fullfile(outputFolder,'Result_RhoOnly.mat');

if ~exist(fileSIMP,'file')
    error('Arquivo nao encontrado: %s',fileSIMP);
end
if ~exist(fileB,'file')
    error('Arquivo nao encontrado: %s',fileB);
end
if ~exist(fileRho,'file')
    error('Arquivo nao encontrado: %s',fileRho);
end

load(fileSIMP,'resultSIMP');
load(fileB,'resultBOnly');
load(fileRho,'resultRhoOnly');

fprintf('\n============================================================\n');
fprintf('CASE 4/4 - TWO VARIABLES\n');
fprintf('============================================================\n');

caseTwoVariable = TutorialFirstTwovariable();

resultTwoVariable.name = 'Two variables';
resultTwoVariable.tutorial = 'TutorialFirstTwovariable';

resultTwoVariable.complianceHistory = ...
    caseTwoVariable.getComplianceHistory();

resultTwoVariable.costHistory = ...
    caseTwoVariable.getCostHistory();

resultTwoVariable.volumeConstraintHistory = ...
    caseTwoVariable.getVolumeConstraintHistory();

fieldsTwo = caseTwoVariable.getDesignVariableValues();

resultTwoVariable.b         = fieldsTwo.b;
resultTwoVariable.rho       = fieldsTwo.rho;
resultTwoVariable.bSmooth   = fieldsTwo.bSmooth;
resultTwoVariable.rhoSmooth = fieldsTwo.rhoSmooth;

resultTwoVariable.mesh = caseTwoVariable.getMeshData();

resultTwoVariable.settings.volumeTarget = 0.4;
resultTwoVariable.settings.bBounds       = [-0.6,0.6];
resultTwoVariable.settings.rhoBounds     = [1e-3,0.95];

save(fullfile(outputFolder,'Result_TwoVariable.mat'), ...
    'resultTwoVariable','-v7.3');

fprintf('Two-variable completed.\n');
fprintf('Iterations              = %d\n', ...
    numel(resultTwoVariable.complianceHistory)-1);
fprintf('Final compliance        = %.12e\n', ...
    resultTwoVariable.complianceHistory(end));
fprintf('Final volume constraint = %.12e\n', ...
    resultTwoVariable.volumeConstraintHistory(end));

%% Historicos

cSIMP = resultSIMP.complianceHistory(:);
cB    = resultBOnly.complianceHistory(:);
cRho  = resultRhoOnly.complianceHistory(:);
cTwo  = resultTwoVariable.complianceHistory(:);

gSIMP = resultSIMP.volumeConstraintHistory(:);
gRho  = resultRhoOnly.volumeConstraintHistory(:);
gTwo  = resultTwoVariable.volumeConstraintHistory(:);

Cref = cSIMP(1);

if isempty(Cref) || ~isfinite(Cref) || Cref <= 0
    error('Cref invalida.');
end

cSIMPNormalized = cSIMP/Cref;
cBNormalized    = cB/Cref;
cRhoNormalized  = cRho/Cref;
cTwoNormalized  = cTwo/Cref;

itSIMP = (0:numel(cSIMP)-1).';
itB    = (0:numel(cB)-1).';
itRho  = (0:numel(cRho)-1).';
itTwo  = (0:numel(cTwo)-1).';

itGSIMP = (0:numel(gSIMP)-1).';
itGRho  = (0:numel(gRho)-1).';
itGTwo  = (0:numel(gTwo)-1).';

%% Salvar historicos combinados

comparison.reference.compliance = Cref;
comparison.reference.description = ...
    'Initial fully solid SIMP design, rho = 1';

comparison.SIMP.compliancePhysical   = cSIMP;
comparison.SIMP.complianceNormalized = cSIMPNormalized;
comparison.SIMP.cost                 = resultSIMP.costHistory(:);
comparison.SIMP.volumeConstraint     = gSIMP;

comparison.BOnly.compliancePhysical   = cB;
comparison.BOnly.complianceNormalized = cBNormalized;
comparison.BOnly.cost                 = resultBOnly.costHistory(:);

comparison.RhoOnly.compliancePhysical   = cRho;
comparison.RhoOnly.complianceNormalized = cRhoNormalized;
comparison.RhoOnly.cost                 = resultRhoOnly.costHistory(:);
comparison.RhoOnly.volumeConstraint     = gRho;

comparison.TwoVariable.compliancePhysical   = cTwo;
comparison.TwoVariable.complianceNormalized = cTwoNormalized;
comparison.TwoVariable.cost = ...
    resultTwoVariable.costHistory(:);
comparison.TwoVariable.volumeConstraint = gTwo;

comparison.settings.volumeTarget  = 0.4;
comparison.settings.bBounds       = [-0.6,0.6];
comparison.settings.normalization = 'Cnormalized = C/Cref';

save(fullfile(outputFolder,'TopologyComparisonHistories.mat'), ...
    'comparison','-v7.3');

%% Grafico 1: compliance normalizada

figCompliance = figure('Color','w');

plot(itSIMP,cSIMPNormalized,'LineWidth',2, ...
    'DisplayName','SIMP');
hold on;
plot(itB,cBNormalized,'LineWidth',2, ...
    'DisplayName','$b$ variable, $\rho=0.4$');
plot(itRho,cRhoNormalized,'LineWidth',2, ...
    'DisplayName','$b=0$, $\rho$ variable');
plot(itTwo,cTwoNormalized,'LineWidth',2, ...
    'DisplayName','$b$ and $\rho$ variable');

yline(1,'--','$C/C_{\mathrm{ref}}=1$', ...
    'Interpreter','latex');

grid on;
box on;
xlabel('Iteration');
ylabel('$C/C_{\mathrm{ref}}$','Interpreter','latex');
title('Normalized compliance comparison');
legend('Interpreter','latex','Location','best');
set(gca,'FontSize',12);

exportgraphics(figCompliance, ...
    fullfile(outputFolder,'NormalizedComplianceComparison.png'), ...
    'Resolution',300);
exportgraphics(figCompliance, ...
    fullfile(outputFolder,'NormalizedComplianceComparison.pdf'), ...
    'ContentType','vector');
savefig(figCompliance, ...
    fullfile(outputFolder,'NormalizedComplianceComparison.fig'));

%% Grafico 2: compliance fisica

figPhysical = figure('Color','w');

plot(itSIMP,cSIMP,'LineWidth',2,'DisplayName','SIMP');
hold on;
plot(itB,cB,'LineWidth',2, ...
    'DisplayName','$b$ variable, $\rho=0.4$');
plot(itRho,cRho,'LineWidth',2, ...
    'DisplayName','$b=0$, $\rho$ variable');
plot(itTwo,cTwo,'LineWidth',2, ...
    'DisplayName','$b$ and $\rho$ variable');

grid on;
box on;
xlabel('Iteration');
ylabel('Compliance');
title('Physical compliance comparison');
legend('Interpreter','latex','Location','best');
set(gca,'FontSize',12);

exportgraphics(figPhysical, ...
    fullfile(outputFolder,'PhysicalComplianceComparison.png'), ...
    'Resolution',300);
exportgraphics(figPhysical, ...
    fullfile(outputFolder,'PhysicalComplianceComparison.pdf'), ...
    'ContentType','vector');
savefig(figPhysical, ...
    fullfile(outputFolder,'PhysicalComplianceComparison.fig'));

%% Grafico 3: volume constraint

figVolume = figure('Color','w');

plot(itGSIMP,gSIMP,'LineWidth',2,'DisplayName','SIMP');
hold on;
plot(itGRho,gRho,'LineWidth',2, ...
    'DisplayName','$b=0$, $\rho$ variable');
plot(itGTwo,gTwo,'LineWidth',2, ...
    'DisplayName','$b$ and $\rho$ variable');

yline(0,'--','Constraint satisfied');

grid on;
box on;
xlabel('Iteration');
ylabel('Volume constraint');
title('Volume-constraint comparison');
legend('Interpreter','latex','Location','best');
set(gca,'FontSize',12);

exportgraphics(figVolume, ...
    fullfile(outputFolder,'VolumeConstraintComparison.png'), ...
    'Resolution',300);
exportgraphics(figVolume, ...
    fullfile(outputFolder,'VolumeConstraintComparison.pdf'), ...
    'ContentType','vector');
savefig(figVolume, ...
    fullfile(outputFolder,'VolumeConstraintComparison.fig'));

%% Tabela final

Case = {
    'SIMP';
    'B-only';
    'Rho-only';
    'Two variables'
};

InitialCompliance = [
    cSIMP(1);
    cB(1);
    cRho(1);
    cTwo(1)
];

FinalCompliance = [
    cSIMP(end);
    cB(end);
    cRho(end);
    cTwo(end)
];

InitialNormalizedCompliance = [
    cSIMPNormalized(1);
    cBNormalized(1);
    cRhoNormalized(1);
    cTwoNormalized(1)
];

FinalNormalizedCompliance = [
    cSIMPNormalized(end);
    cBNormalized(end);
    cRhoNormalized(end);
    cTwoNormalized(end)
];

FinalVolumeConstraint = [
    gSIMP(end);
    NaN;
    gRho(end);
    gTwo(end)
];

Iterations = [
    numel(cSIMP)-1;
    numel(cB)-1;
    numel(cRho)-1;
    numel(cTwo)-1
];

summaryTable = table( ...
    Case, ...
    InitialCompliance, ...
    FinalCompliance, ...
    InitialNormalizedCompliance, ...
    FinalNormalizedCompliance, ...
    FinalVolumeConstraint, ...
    Iterations);

disp(summaryTable);

writetable(summaryTable, ...
    fullfile(outputFolder,'TopologyComparisonSummary.csv'));

save(fullfile(outputFolder,'TopologyComparisonSummary.mat'), ...
    'summaryTable','Cref');

fprintf('\n============================================================\n');
fprintf('LAST CASE AND COMPARISON COMPLETED\n');
fprintf('Results saved in:\n%s\n',outputFolder);
fprintf('============================================================\n');
