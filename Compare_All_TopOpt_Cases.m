clear
clc
close all

outFolder = 'Comparison_TopOpt_4Cases';

if ~exist(outFolder,'dir')
    mkdir(outFolder);
end

A = load('Result_A_Rho.mat');
B = load('Result_B_Rho.mat');
C = load('Result_RhoOnly.mat');
D = load('Result_SIMP.mat');

ARho    = A.R;
BRho    = B.R;
RhoOnly = C.R;
SIMP    = D.R;

J1 = ARho.compliance(:);
J2 = BRho.compliance(:);
J3 = RhoOnly.compliance(:);
J4 = SIMP.compliance(:);

it1 = 0:length(J1)-1;
it2 = 0:length(J2)-1;
it3 = 0:length(J3)-1;
it4 = 0:length(J4)-1;

fig = figure( ...
    'Color','w', ...
    'Position',[100 100 900 650]);

plot(it1,J1/J1(1),'LineWidth',2);
hold on
plot(it2,J2/J2(1),'LineWidth',2);
plot(it3,J3/J3(1),'LineWidth',2);
plot(it4,J4/J4(1),'LineWidth',2);

xlabel('Iteration');
ylabel('Normalized compliance');

legend( ...
    'a + \rho', ...
    'b + \rho', ...
    '\rho only', ...
    'SIMP', ...
    'Location','best');

grid on
box on

exportgraphics( ...
    fig, ...
    fullfile(outFolder,'Compliance_Comparison_4Cases.png'), ...
    'Resolution',300);

close(fig)

g1 = ARho.volume(:);
g2 = BRho.volume(:);
g3 = RhoOnly.volume(:);
g4 = SIMP.volume(:);

it1 = 0:length(g1)-1;
it2 = 0:length(g2)-1;
it3 = 0:length(g3)-1;
it4 = 0:length(g4)-1;

fig = figure( ...
    'Color','w', ...
    'Position',[100 100 900 650]);

plot(it1,g1,'LineWidth',2);
hold on
plot(it2,g2,'LineWidth',2);
plot(it3,g3,'LineWidth',2);
plot(it4,g4,'LineWidth',2);

yline(0,'k--','LineWidth',1.5);

xlabel('Iteration');
ylabel('Volume constraint');

legend( ...
    'a + \rho', ...
    'b + \rho', ...
    '\rho only', ...
    'SIMP', ...
    'Target', ...
    'Location','best');

grid on
box on

exportgraphics( ...
    fig, ...
    fullfile(outFolder,'Volume_Comparison_4Cases.png'), ...
    'Resolution',300);

close(fig)

plotField( ...
    ARho.mesh, ...
    ARho.rhoSmooth, ...
    '\rho - a+\rho', ...
    fullfile(outFolder,'Field_Rho_A_Rho.png'));

plotField( ...
    ARho.mesh, ...
    ARho.aSmooth, ...
    'a - a+\rho', ...
    fullfile(outFolder,'Field_A_A_Rho.png'));

plotField( ...
    BRho.mesh, ...
    BRho.rhoSmooth, ...
    '\rho - b+\rho', ...
    fullfile(outFolder,'Field_Rho_B_Rho.png'));

plotField( ...
    BRho.mesh, ...
    BRho.bSmooth, ...
    'b - b+\rho', ...
    fullfile(outFolder,'Field_B_B_Rho.png'));

plotField( ...
    RhoOnly.mesh, ...
    RhoOnly.rhoSmooth, ...
    '\rho - rho only', ...
    fullfile(outFolder,'Field_RhoOnly.png'));

plotField( ...
    SIMP.mesh, ...
    SIMP.rhoSmooth, ...
    '\rho - SIMP', ...
    fullfile(outFolder,'Field_SIMP.png'));

fprintf('\n============================================\n');
fprintf('FINAL RESULTS\n');
fprintf('============================================\n');

fprintf('\na + rho:\n');
fprintf('  compliance = %.12e\n',J1(end));
fprintf('  volume     = %.12e\n',g1(end));

fprintf('\nb + rho:\n');
fprintf('  compliance = %.12e\n',J2(end));
fprintf('  volume     = %.12e\n',g2(end));

fprintf('\nrho only:\n');
fprintf('  compliance = %.12e\n',J3(end));
fprintf('  volume     = %.12e\n',g3(end));

fprintf('\nSIMP:\n');
fprintf('  compliance = %.12e\n',J4(end));
fprintf('  volume     = %.12e\n',g4(end));

fprintf('\n============================================\n');
fprintf('B FIELD\n');
fprintf('============================================\n');

fprintf('b min       = %.8f\n',min(BRho.b));
fprintf('b max       = %.8f\n',max(BRho.b));
fprintf('bSmooth min = %.8f\n',min(BRho.bSmooth));
fprintf('bSmooth max = %.8f\n',max(BRho.bSmooth));

fprintf('\nResults saved in:\n');
fprintf('%s\n',outFolder);

function plotField(meshData,values,ttl,fileName)

    values = values(:);

    fig = figure( ...
        'Color','w', ...
        'Position',[100 100 1000 500]);

    patch( ...
        'Faces',meshData.connec, ...
        'Vertices',meshData.coord, ...
        'FaceVertexCData',values, ...
        'FaceColor','interp', ...
        'EdgeColor','none');

    axis equal
    axis tight

    xlabel('x');
    ylabel('y');

    title(ttl);

    colorbar
    box on

    exportgraphics( ...
        fig, ...
        fileName, ...
        'Resolution',300);

    close(fig)

end