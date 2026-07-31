%% RunTopologyComparison
% Executa quatro casos de topology optimization:
% 1) SIMP
% 2) B-only
% 3) Rho-only
% 4) B e rho variaveis
%
% Referencia comum:
% Cref = compliance inicial do caso SIMP totalmente solido.
%
% Observacao:
% O caso B-only nao possui volume constraint.

clear;
close all;
clc;

clear classes;
rehash;

%% Pasta de resultados

outputFolder = fullfile(pwd,'TopologyComparisonResults');

if ~exist(outputFolder,'dir')
    mkdir(outputFolder);
end

fprintf('\n');
fprintf('============================================================\n');
fprintf('TOPOLOGY OPTIMIZATION COMPARISON\n');
fprintf('============================================================\n');
fprintf('Output folder:\n%s\n',outputFolder);
fprintf('============================================================\n');


%% ============================================================
%  CASO 1: SIMP
%  ============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('CASE 1/4 - SIMP REFERENCE\n');
fprintf('Tutorial05_2_TopOpt2DDensityMacroNullSpace\n');
fprintf('============================================================\n');

caseSIMP = Tutorial05_2_TopOpt2DDensityMacroNullSpace();

resultSIMP.name = 'SIMP';
resultSIMP.tutorial = ...
    'Tutorial05_2_TopOpt2DDensityMacroNullSpace';

resultSIMP.complianceHistory = ...
    caseSIMP.getComplianceHistory();

resultSIMP.costHistory = ...
    caseSIMP.getCostHistory();

resultSIMP.volumeConstraintHistory = ...
    caseSIMP.getVolumeConstraintHistory();

% Campo rho final
resultSIMP.rho = ...
    caseSIMP.getRhoValues();

% Campo rho filtrado final
resultSIMP.rhoSmooth = ...
    caseSIMP.getRhoSmoothValues();

% Dados numericos da malha
resultSIMP.mesh = ...
    caseSIMP.getMeshData();

resultSIMP.settings.volumeTarget = 0.4;
resultSIMP.settings.rhoBounds    = [0,1];

save( ...
    fullfile(outputFolder,'Result_SIMP.mat'), ...
    'resultSIMP', ...
    '-v7.3');

fprintf('SIMP completed.\n');
fprintf('Iterations              = %d\n', ...
    numel(resultSIMP.complianceHistory)-1);
fprintf('Final compliance        = %.12e\n', ...
    resultSIMP.complianceHistory(end));
fprintf('Final volume constraint = %.12e\n', ...
    resultSIMP.volumeConstraintHistory(end));


%% ============================================================
%  CASO 2: B-ONLY
%  ============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('CASE 2/4 - B-ONLY\n');
fprintf('TutorialFirstBOnly\n');
fprintf('============================================================\n');

caseBOnly = TutorialFirstBOnly();

resultBOnly.name = 'B-only';
resultBOnly.tutorial = 'TutorialFirstBOnly';

resultBOnly.complianceHistory = ...
    caseBOnly.getComplianceHistory();

resultBOnly.costHistory = ...
    caseBOnly.getCostHistory();

resultBOnly.beta = ...
    caseBOnly.getBValues();

resultBOnly.betaSmooth = ...
    caseBOnly.getBSmoothValues();

resultBOnly.bGeomSmooth = ...
    caseBOnly.getBGeomSmoothValues();

resultBOnly.mesh = ...
    caseBOnly.getMeshData();

resultBOnly.settings.fixedRho    = 0.4;
resultBOnly.settings.betaBounds  = [-0.6,0.6];
resultBOnly.settings.bGeomBounds = [0,0.6];

save( ...
    fullfile(outputFolder,'Result_BOnly.mat'), ...
    'resultBOnly', ...
    '-v7.3');

fprintf('B-only completed.\n');
fprintf('Iterations       = %d\n', ...
    numel(resultBOnly.complianceHistory)-1);
fprintf('Final compliance = %.12e\n', ...
    resultBOnly.complianceHistory(end));


%% ============================================================
%  CASO 3: RHO-ONLY
%  ============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('CASE 3/4 - RHO-ONLY\n');
fprintf('TutorialFirstRhoOnly\n');
fprintf('============================================================\n');

caseRhoOnly = TutorialFirstRhoOnly();

resultRhoOnly.name = 'Rho-only';
resultRhoOnly.tutorial = 'TutorialFirstRhoOnly';

resultRhoOnly.complianceHistory = ...
    caseRhoOnly.getComplianceHistory();

resultRhoOnly.costHistory = ...
    caseRhoOnly.getCostHistory();

resultRhoOnly.volumeConstraintHistory = ...
    caseRhoOnly.getVolumeConstraintHistory();

resultRhoOnly.rho = ...
    caseRhoOnly.getRhoValues();

resultRhoOnly.rhoSmooth = ...
    caseRhoOnly.getRhoSmoothValues();

resultRhoOnly.b = ...
    caseRhoOnly.getBValues();

resultRhoOnly.mesh = ...
    caseRhoOnly.getMeshData();

resultRhoOnly.settings.fixedB       = 0;
resultRhoOnly.settings.volumeTarget = 0.4;
resultRhoOnly.settings.rhoBounds    = [1e-3,0.97];

save( ...
    fullfile(outputFolder,'Result_RhoOnly.mat'), ...
    'resultRhoOnly', ...
    '-v7.3');

fprintf('Rho-only completed.\n');
fprintf('Iterations              = %d\n', ...
    numel(resultRhoOnly.complianceHistory)-1);
fprintf('Final compliance        = %.12e\n', ...
    resultRhoOnly.complianceHistory(end));
fprintf('Final volume constraint = %.12e\n', ...
    resultRhoOnly.volumeConstraintHistory(end));


%% ============================================================
%  CASO 4: TWO VARIABLES
%  ============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('CASE 4/4 - TWO VARIABLES\n');
fprintf('TutorialFirstTwovariable\n');
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

fieldsTwo = ...
    caseTwoVariable.getDesignVariableValues();

resultTwoVariable.b = fieldsTwo.b;
resultTwoVariable.rho = fieldsTwo.rho;
resultTwoVariable.bSmooth = fieldsTwo.bSmooth;
resultTwoVariable.rhoSmooth = fieldsTwo.rhoSmooth;

resultTwoVariable.mesh = ...
    caseTwoVariable.getMeshData();

resultTwoVariable.settings.volumeTarget = 0.4;
resultTwoVariable.settings.bBounds       = [-0.6,0.6];
resultTwoVariable.settings.rhoBounds     = [1e-3,0.95];

save( ...
    fullfile(outputFolder,'Result_TwoVariable.mat'), ...
    'resultTwoVariable', ...
    '-v7.3');

fprintf('Two-variable case completed.\n');
fprintf('Iterations              = %d\n', ...
    numel(resultTwoVariable.complianceHistory)-1);
fprintf('Final compliance        = %.12e\n', ...
    resultTwoVariable.complianceHistory(end));
fprintf('Final volume constraint = %.12e\n', ...
    resultTwoVariable.volumeConstraintHistory(end));


%% ============================================================
%  EXTRAIR HISTORICOS
%  ============================================================

cSIMP = resultSIMP.complianceHistory(:);
cB    = resultBOnly.complianceHistory(:);
cRho  = resultRhoOnly.complianceHistory(:);
cTwo  = resultTwoVariable.complianceHistory(:);

gSIMP = resultSIMP.volumeConstraintHistory(:);
gRho  = resultRhoOnly.volumeConstraintHistory(:);
gTwo  = resultTwoVariable.volumeConstraintHistory(:);


%% ============================================================
%  COMPLIANCE DE REFERENCIA
%  ============================================================

% O SIMP inicia com rho = 1 em todo o dominio.
Cref = cSIMP(1);

if isempty(Cref)
    error('A compliance de referencia esta vazia.');
end

if ~isfinite(Cref)
    error('A compliance de referencia nao e finita.');
end

if Cref <= 0
    error('A compliance de referencia deve ser positiva.');
end

fprintf('\n');
fprintf('============================================================\n');
fprintf('COMMON COMPLIANCE REFERENCE\n');
fprintf('============================================================\n');
fprintf('Cref = initial SIMP compliance\n');
fprintf('Cref = %.12e\n',Cref);
fprintf('============================================================\n');


%% ============================================================
%  NORMALIZAR COMPLIANCE
%  ============================================================

cSIMPNormalized = cSIMP/Cref;
cBNormalized    = cB/Cref;
cRhoNormalized  = cRho/Cref;
cTwoNormalized  = cTwo/Cref;


%% ============================================================
%  VETORES DE ITERACAO
%  ============================================================

itSIMP = (0:numel(cSIMP)-1).';
itB    = (0:numel(cB)-1).';
itRho  = (0:numel(cRho)-1).';
itTwo  = (0:numel(cTwo)-1).';

itGSIMP = (0:numel(gSIMP)-1).';
itGRho  = (0:numel(gRho)-1).';
itGTwo  = (0:numel(gTwo)-1).';


%% ============================================================
%  STRUCT COM TODOS OS HISTORICOS
%  ============================================================

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
comparison.TwoVariable.cost                 = resultTwoVariable.costHistory(:);
comparison.TwoVariable.volumeConstraint     = gTwo;

comparison.settings.volumeTarget  = 0.4;
comparison.settings.bBounds       = [-0.6,0.6];
comparison.settings.normalization = 'Cnormalized = C/Cref';

save( ...
    fullfile(outputFolder,'TopologyComparisonHistories.mat'), ...
    'comparison', ...
    '-v7.3');


%% ============================================================
%  GRAFICO 1: COMPLIANCE NORMALIZADA
%  ============================================================

figCompliance = figure( ...
    'Color','w', ...
    'Name','Normalized compliance comparison');

plot(itSIMP,cSIMPNormalized, ...
    'LineWidth',2, ...
    'DisplayName','SIMP');

hold on;

plot(itB,cBNormalized, ...
    'LineWidth',2, ...
    'DisplayName','$b$ variable, $\rho=0.4$');

plot(itRho,cRhoNormalized, ...
    'LineWidth',2, ...
    'DisplayName','$b=0$, $\rho$ variable');

plot(itTwo,cTwoNormalized, ...
    'LineWidth',2, ...
    'DisplayName','$b$ and $\rho$ variable');

yline(1,'--', ...
    '$C/C_{\mathrm{ref}}=1$', ...
    'Interpreter','latex', ...
    'LabelHorizontalAlignment','left');

grid on;
box on;

xlabel('Iteration');
ylabel('$C/C_{\mathrm{ref}}$','Interpreter','latex');
title('Normalized compliance comparison');

legend( ...
    'Interpreter','latex', ...
    'Location','best');

set(gca,'FontSize',12);

exportgraphics( ...
    figCompliance, ...
    fullfile(outputFolder,'NormalizedComplianceComparison.png'), ...
    'Resolution',300);

exportgraphics( ...
    figCompliance, ...
    fullfile(outputFolder,'NormalizedComplianceComparison.pdf'), ...
    'ContentType','vector');

savefig( ...
    figCompliance, ...
    fullfile(outputFolder,'NormalizedComplianceComparison.fig'));


%% ============================================================
%  GRAFICO EXTRA: COMPLIANCE FISICA
%  ============================================================

figCompliancePhysical = figure( ...
    'Color','w', ...
    'Name','Physical compliance comparison');

plot(itSIMP,cSIMP, ...
    'LineWidth',2, ...
    'DisplayName','SIMP');

hold on;

plot(itB,cB, ...
    'LineWidth',2, ...
    'DisplayName','$b$ variable, $\rho=0.4$');

plot(itRho,cRho, ...
    'LineWidth',2, ...
    'DisplayName','$b=0$, $\rho$ variable');

plot(itTwo,cTwo, ...
    'LineWidth',2, ...
    'DisplayName','$b$ and $\rho$ variable');

grid on;
box on;

xlabel('Iteration');
ylabel('Compliance');
title('Physical compliance comparison');

legend( ...
    'Interpreter','latex', ...
    'Location','best');

set(gca,'FontSize',12);

exportgraphics( ...
    figCompliancePhysical, ...
    fullfile(outputFolder,'PhysicalComplianceComparison.png'), ...
    'Resolution',300);

exportgraphics( ...
    figCompliancePhysical, ...
    fullfile(outputFolder,'PhysicalComplianceComparison.pdf'), ...
    'ContentType','vector');

savefig( ...
    figCompliancePhysical, ...
    fullfile(outputFolder,'PhysicalComplianceComparison.fig'));


%% ============================================================
%  GRAFICO 2: VOLUME CONSTRAINT
%  ============================================================

% B-only nao aparece porque rho = 0.4 e fixo.

figVolume = figure( ...
    'Color','w', ...
    'Name','Volume constraint comparison');

plot(itGSIMP,gSIMP, ...
    'LineWidth',2, ...
    'DisplayName','SIMP');

hold on;

plot(itGRho,gRho, ...
    'LineWidth',2, ...
    'DisplayName','$b=0$, $\rho$ variable');

plot(itGTwo,gTwo, ...
    'LineWidth',2, ...
    'DisplayName','$b$ and $\rho$ variable');

yline(0,'--', ...
    'Constraint satisfied', ...
    'LabelHorizontalAlignment','left');

grid on;
box on;

xlabel('Iteration');
ylabel('Volume constraint');
title('Volume-constraint comparison');

legend( ...
    'Interpreter','latex', ...
    'Location','best');

set(gca,'FontSize',12);

exportgraphics( ...
    figVolume, ...
    fullfile(outputFolder,'VolumeConstraintComparison.png'), ...
    'Resolution',300);

exportgraphics( ...
    figVolume, ...
    fullfile(outputFolder,'VolumeConstraintComparison.pdf'), ...
    'ContentType','vector');

savefig( ...
    figVolume, ...
    fullfile(outputFolder,'VolumeConstraintComparison.fig'));


%% ============================================================
%  TABELA FINAL
%  ============================================================

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

disp(' ');
disp('============================================================');
disp('FINAL COMPARISON');
disp('============================================================');
disp(summaryTable);

fprintf('\n');
fprintf('Common compliance reference:\n');
fprintf('Cref = %.12e\n',Cref);

writetable( ...
    summaryTable, ...
    fullfile(outputFolder,'TopologyComparisonSummary.csv'));

save( ...
    fullfile(outputFolder,'TopologyComparisonSummary.mat'), ...
    'summaryTable', ...
    'Cref');


%% ============================================================
%  FINALIZACAO
%  ============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('ALL CASES COMPLETED\n');
fprintf('============================================================\n');
fprintf('Results saved in:\n');
fprintf('%s\n',outputFolder);
fprintf('============================================================\n');
