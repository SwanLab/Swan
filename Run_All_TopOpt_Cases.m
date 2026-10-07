clear
clc
close all

outFolder = 'Comparison_TopOpt';

if ~exist(outFolder,'dir')
    mkdir(outFolder);
end


fprintf('\n============================================\n');
fprintf('CASE 1: a + rho\n');
fprintf('============================================\n');

objARho = TutorialFirstTwovariable();

R = struct();

R.caseName   = 'A_Rho';

R.cost       = objARho.getCostHistory();
R.compliance = objARho.getComplianceHistory();
R.volume     = objARho.getVolumeConstraintHistory();

R.a          = objARho.getAValues();
R.aSmooth    = objARho.getASmoothValues();

R.rho        = objARho.getRhoValues();
R.rhoSmooth  = objARho.getRhoSmoothValues();

R.mesh       = objARho.getMeshData();

save('Result_A_Rho.mat','R','-v7.3');

clear objARho R


fprintf('\n============================================\n');
fprintf('CASE 2: rho only\n');
fprintf('============================================\n');

objRho = TutorialFirstRhoOnly();

R = struct();

R.caseName   = 'RhoOnly';

R.cost       = objRho.getCostHistory();
R.compliance = objRho.getComplianceHistory();
R.volume     = objRho.getVolumeConstraintHistory();

R.a          = objRho.getAValues();
R.aSmooth    = objRho.getASmoothValues();

R.rho        = objRho.getRhoValues();
R.rhoSmooth  = objRho.getRhoSmoothValues();

R.mesh       = objRho.getMeshData();

save('Result_RhoOnly.mat','R','-v7.3');

clear objRho R


fprintf('\n============================================\n');
fprintf('CASE 3: SIMP\n');
fprintf('============================================\n');

objSIMP = Tutorial05_2_TopOpt2DDensityMacroNullSpace();

R = struct();

R.caseName   = 'SIMP';

R.cost       = objSIMP.getCostHistory();
R.compliance = objSIMP.getComplianceHistory();
R.volume     = objSIMP.getVolumeConstraintHistory();

R.rho        = objSIMP.getRhoValues();
R.rhoSmooth  = objSIMP.getRhoSmoothValues();

R.mesh       = objSIMP.getMeshData();

save('Result_SIMP.mat','R','-v7.3');

clear objSIMP R


fprintf('\n============================================\n');
fprintf('POST-PROCESSING\n');
fprintf('============================================\n');

A = load('Result_A_Rho.mat');
B = load('Result_RhoOnly.mat');
C = load('Result_SIMP.mat');

ARho    = A.R;
RhoOnly = B.R;
SIMP    = C.R;


J1 = ARho.compliance(:);
J2 = RhoOnly.compliance(:);
J3 = SIMP.compliance(:);

it1 = 0:length(J1)-1;
it2 = 0:length(J2)-1;
it3 = 0:length(J3)-1;


fig = figure( ...
    'Color','w', ...
    'Position',[100 100 900 650]);

plot(it1,J1/J1(1),'LineWidth',2);
hold on

plot(it2,J2/J2(1),'LineWidth',2);
plot(it3,J3/J3(1),'LineWidth',2);

xlabel('Iteration');
ylabel('Normalized compliance');

legend( ...
    'a + \rho', ...
    '\rho only', ...
    'SIMP', ...
    'Location','best');

grid on
box on

exportgraphics( ...
    fig, ...
    fullfile(outFolder,'Compliance_Comparison.png'), ...
    'Resolution',300);

close(fig)


g1 = ARho.volume(:);
g2 = RhoOnly.volume(:);
g3 = SIMP.volume(:);

it1 = 0:length(g1)-1;
it2 = 0:length(g2)-1;
it3 = 0:length(g3)-1;


fig = figure( ...
    'Color','w', ...
    'Position',[100 100 900 650]);

plot(it1,g1,'LineWidth',2);
hold on

plot(it2,g2,'LineWidth',2);
plot(it3,g3,'LineWidth',2);

yline(0,'k--','LineWidth',1.5);

xlabel('Iteration');
ylabel('Volume constraint');

legend( ...
    'a + \rho', ...
    '\rho only', ...
    'SIMP', ...
    'Target', ...
    'Location','best');

grid on
box on

exportgraphics( ...
    fig, ...
    fullfile(outFolder,'Volume_Comparison.png'), ...
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

fprintf('a + rho:\n');
fprintf('  compliance = %.12e\n',J1(end));
fprintf('  volume     = %.12e\n',g1(end));

fprintf('\nrho only:\n');
fprintf('  compliance = %.12e\n',J2(end));
fprintf('  volume     = %.12e\n',g2(end));

fprintf('\nSIMP:\n');
fprintf('  compliance = %.12e\n',J3(end));
fprintf('  volume     = %.12e\n',g3(end));

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