clear; close all; clc;
outDir = 'C:\Users\david\Desktop\PlotsResultats\1.Single\ModeI';

if ~exist(outDir,'dir')
    mkdir(outDir);
end

figPos = [100 100 900 650];
fontSize = 18;
lineWidth = 1.5;

%% Load data

dataTotal = load('SingleTotal15_30.mat');
dataUnload = load('SingleUnload15_30.mat');

Ftot = dataTotal.Res(:,1);
utot = dataTotal.Res(:,2)*1e3;

Funl = dataUnload.Res(:,1);
uunl = dataUnload.Res(:,2)*1e3;

dtot = dataTotal.Res(:,3);
dunl = dataUnload.Res(:,3);

idxComp = find(uunl < 0);
Funl(idxComp) = -Funl(idxComp);

%% Compute external work (area under force-displacement curve)

Wtot = cumtrapz(utot,Ftot)*1e-3;

%% Force-displacement (total loading)

figure;
plot(utot,Ftot,'LineWidth',1.5);
set(gca,'FontSize',fontSize)
grid on;
xlabel('Displacement [mm]');
ylabel('Force [N]');
title('Single Element Benchmark - Total Loading');

%% Force-displacement (total loading)

fig = figure('Position',figPos);
plot(utot,Ftot,'LineWidth',1.5);
set(gca,'FontSize',fontSize)
grid on;
xlabel('Displacement [mm]');
ylabel('Force [N]');
title('Single Element Benchmark - Total Loading');

exportgraphics(fig,fullfile(outDir,'ForceDisplacement_Total.png'),'Resolution',300);

%% Force-displacement (loading/unloading)

fig = figure('Position',figPos);
plot(utot,Ftot,'LineWidth',1.5);
hold on;
plot(uunl,Funl,'LineWidth',1.5);
set(gca,'FontSize',fontSize)
grid on;
xlabel('Displacement [mm]');
ylabel('Force [N]');
title('Single Element Benchmark - Loading and Unloading');
legend('Loading','Unloading','Location','best');

exportgraphics(fig,fullfile(outDir,'ForceDisplacement_LoadUnload.png'),'Resolution',300);

%% Energy evolution

fig = figure('Position',figPos);
plot(utot,Wtot,'LineWidth',1.5);
set(gca,'FontSize',fontSize)
grid on;
xlabel('Displacement [mm]');
ylabel('Energy [J]');
title('Accumulated External Work');

exportgraphics(fig,fullfile(outDir,'EnergyEvolution.png'),'Resolution',300);

%% Damage evolution (total loading)

fig = figure('Position',figPos);
plot(utot,dtot,'LineWidth',1.5);
set(gca,'FontSize',fontSize)
grid on;
xlabel('Displacement');
ylabel('Damage');
title('Damage Evolution - Total Loading');

exportgraphics(fig,fullfile(outDir,'Damage_Total.png'),'Resolution',300);

%% Damage evolution (loading/unloading)

fig = figure('Position',figPos);
plot(utot,dtot,'LineWidth',1.5);
set(gca,'FontSize',fontSize)
hold on;
plot(uunl,dunl,'LineWidth',1.5);
set(gca,'FontSize',fontSize)
grid on;
xlabel('Displacement');
ylabel('Damage');
title('Damage Evolution - Loading and Unloading');
legend('Loading','Unloading','Location','best');

exportgraphics(fig,fullfile(outDir,'Damage_LoadUnload.png'),'Resolution',300);

%% Print final energy

fprintf('Final energy = %.6e\n',Wtot(end));