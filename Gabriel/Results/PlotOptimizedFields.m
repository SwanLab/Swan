%% PlotOptimizedFields
clear; close all; clc;

resultsFolder = fullfile(pwd,'TopologyComparisonResults');

load(fullfile(resultsFolder,'Result_SIMP.mat'),'resultSIMP');
load(fullfile(resultsFolder,'Result_BOnly.mat'),'resultBOnly');
load(fullfile(resultsFolder,'Result_RhoOnly.mat'),'resultRhoOnly');
load(fullfile(resultsFolder,'Result_TwoVariable.mat'),'resultTwoVariable');

plotNodalField(resultSIMP.mesh,resultSIMP.rhoSmooth, ...
    'SIMP: optimized density field $\\rho$',[0 1]);
saveCurrent(resultsFolder,'Field_SIMP_RhoSmooth');

plotNodalField(resultBOnly.mesh,resultBOnly.betaSmooth, ...
    'B-only: optimized field $\\beta$',[-0.6 0.6]);
saveCurrent(resultsFolder,'Field_BOnly_BetaSmooth');

plotNodalField(resultBOnly.mesh,resultBOnly.bGeomSmooth, ...
    'B-only: geometric parameter $b_{\\mathrm{geom}}$',[0 0.6]);
saveCurrent(resultsFolder,'Field_BOnly_BGeomSmooth');

plotNodalField(resultRhoOnly.mesh,resultRhoOnly.rhoSmooth, ...
    'Rho-only: optimized density field $\\rho$',[0 1]);
saveCurrent(resultsFolder,'Field_RhoOnly_RhoSmooth');

plotNodalField(resultTwoVariable.mesh,resultTwoVariable.bSmooth, ...
    'Two variables: optimized field $b$',[-0.6 0.6]);
saveCurrent(resultsFolder,'Field_TwoVariable_BSmooth');

plotNodalField(resultTwoVariable.mesh,resultTwoVariable.rhoSmooth, ...
    'Two variables: optimized density field $\\rho$',[0 1]);
saveCurrent(resultsFolder,'Field_TwoVariable_RhoSmooth');

fprintf('Optimized fields plotted and saved in:\n%s\n',resultsFolder);

function plotNodalField(meshData,fieldData,plotTitle,colorLimits)
coord = meshData.coord;
connec = meshData.connec;
if size(connec,2) > 3
    connec = connec(:,1:3);
end
fieldData = fieldData(:);
if numel(fieldData) ~= size(coord,1)
    error('Field has %d values, mesh has %d nodes.', ...
        numel(fieldData),size(coord,1));
end
figure('Color','w');
patch('Faces',connec,'Vertices',coord(:,1:2), ...
    'FaceVertexCData',fieldData, ...
    'FaceColor','interp','EdgeColor','none');
axis equal tight; box on;
xlabel('$x$','Interpreter','latex');
ylabel('$y$','Interpreter','latex');
title(plotTitle,'Interpreter','latex');
colorbar;
clim(colorLimits);
set(gca,'FontSize',12,'TickLabelInterpreter','latex');
end

function saveCurrent(resultsFolder,fileName)
fig = gcf;
exportgraphics(fig,fullfile(resultsFolder,[fileName '.png']), ...
    'Resolution',300);
savefig(fig,fullfile(resultsFolder,[fileName '.fig']));
end
