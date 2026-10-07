close all
clc

outFolder = 'PostProcess_B_Rho_Output';

if ~exist(outFolder,'dir')
    mkdir(outFolder);
end

H = load('HomogenizationtwovariablesVB.mat');

Chomog        = H.Chomog;
paramB        = H.paramB(:);
paramRho      = H.paramRho(:);
volFrac       = H.volFrac;
Interpolation = H.Interpolation;

nB   = numel(paramB);
nRho = numel(paramRho);

[Bgrid,RHOgrid] = meshgrid(paramB,paramRho);

compNames = { ...
    'C1111', ...
    'C1122', ...
    'C1112', ...
    'C2222', ...
    'C2212', ...
    'C1212'};

compIds = { ...
    [1 1 1 1]
    [1 1 2 2]
    [1 1 1 2]
    [2 2 2 2]
    [2 2 1 2]
    [1 2 1 2]};

Yhom = zeros(nRho,nB,6);

for ic = 1:6

    id = compIds{ic};

    Yhom(:,:,ic) = squeeze( ...
        Chomog(id(1),id(2),id(3),id(4),:,:) );

end

fprintf('\n============================================\n');
fprintf('POST-PROCESS B-RHO\n');
fprintf('============================================\n');
fprintf('Output folder: %s\n',outFolder);
fprintf('Grid size     : %d rho x %d b\n',nRho,nB);

fprintf('\nEvaluating NN surfaces...\n');

Ynn = zeros(nRho,nB,6);

for ic = 1:6

    id = compIds{ic};

    funNN = Interpolation.fun{ ...
        id(1),id(2),id(3),id(4)};

    for ir = 1:nRho

        for ib = 1:nB

            Ynn(ir,ib,ic) = ...
                funNN(paramB(ib),paramRho(ir));

        end

    end

end

fprintf('NN surfaces evaluated.\n');

for ic = 1:6

    fig = figure( ...
        'Color','w', ...
        'Position',[80 80 900 700]);

    surf( ...
        Bgrid, ...
        RHOgrid, ...
        Yhom(:,:,ic), ...
        'EdgeColor','none');

    xlabel('b');
    ylabel('\rho');
    zlabel(compNames{ic});

    title(['HOM - ' compNames{ic}]);

    view(45,30);
    grid on
    box on
    colorbar

    exportgraphics( ...
        fig, ...
        fullfile(outFolder, ...
        ['HOM_Surface_' compNames{ic} '.png']), ...
        'Resolution',300);

    close(fig)

end

for ic = 1:6

    fig = figure( ...
        'Color','w', ...
        'Position',[80 80 900 700]);

    surf( ...
        Bgrid, ...
        RHOgrid, ...
        Ynn(:,:,ic), ...
        'EdgeColor','none');

    xlabel('b');
    ylabel('\rho');
    zlabel(compNames{ic});

    title(['NN - ' compNames{ic}]);

    view(45,30);
    grid on
    box on
    colorbar

    exportgraphics( ...
        fig, ...
        fullfile(outFolder, ...
        ['NN_Surface_' compNames{ic} '.png']), ...
        'Resolution',300);

    close(fig)

end

for ic = 1:6

    Zhom = Yhom(:,:,ic);
    Znn  = Ynn(:,:,ic);

    errAbs = abs(Znn-Zhom);

    scale = max(abs(Zhom),[],'all');

    if scale < 1e-14
        scale = 1;
    end

    errGlobal = errAbs/scale;

    fig = figure( ...
        'Color','w', ...
        'Position',[50 50 1400 900]);

    t = tiledlayout(2,2, ...
        'Padding','compact', ...
        'TileSpacing','compact');

    nexttile

    surf( ...
        Bgrid,RHOgrid,Zhom, ...
        'EdgeColor','none');

    xlabel('b');
    ylabel('\rho');
    zlabel(compNames{ic});

    title('HOM');
    view(45,30);
    grid on
    box on
    colorbar

    nexttile

    surf( ...
        Bgrid,RHOgrid,Znn, ...
        'EdgeColor','none');

    xlabel('b');
    ylabel('\rho');
    zlabel(compNames{ic});

    title('NN');
    view(45,30);
    grid on
    box on
    colorbar

    nexttile

    imagesc( ...
        paramB, ...
        paramRho, ...
        errAbs);

    axis xy

    xlabel('b');
    ylabel('\rho');

    title('|NN - HOM|');

    colorbar

    nexttile

    imagesc( ...
        paramB, ...
        paramRho, ...
        errGlobal);

    axis xy

    xlabel('b');
    ylabel('\rho');

    title('Global relative error');

    colorbar

    title(t,['HOM vs NN - ' compNames{ic}]);

    exportgraphics( ...
        fig, ...
        fullfile(outFolder, ...
        ['HOM_vs_NN_' compNames{ic} '.png']), ...
        'Resolution',300);

    close(fig)

end

rhoTarget = 0.5;

[~,irFix] = min(abs(paramRho-rhoTarget));

rhoFix = paramRho(irFix);

fprintf('\nFixed rho for curves = %.8f\n',rhoFix);

fig = figure( ...
    'Color','w', ...
    'Position',[50 50 1400 1000]);

t = tiledlayout(3,2, ...
    'Padding','compact', ...
    'TileSpacing','compact');

for ic = 1:6

    nexttile

    plot( ...
        paramB, ...
        Yhom(irFix,:,ic), ...
        'k-', ...
        'LineWidth',2);

    hold on

    plot( ...
        paramB, ...
        Ynn(irFix,:,ic), ...
        'r--', ...
        'LineWidth',2);

    xlabel('b');
    ylabel(compNames{ic});

    title(compNames{ic});

    legend( ...
        'HOM', ...
        'NN', ...
        'Location','best');

    grid on
    box on

end

title(t, ...
    sprintf('Tensor components vs b, \\rho = %.4f',rhoFix));

exportgraphics( ...
    fig, ...
    fullfile(outFolder, ...
    'Components_vs_b_rho_050.png'), ...
    'Resolution',300);

close(fig)

for ic = 1:6

    fig = figure( ...
        'Color','w', ...
        'Position',[100 100 850 650]);

    plot( ...
        paramB, ...
        Yhom(irFix,:,ic), ...
        'ko-', ...
        'LineWidth',1.5, ...
        'MarkerSize',4);

    hold on

    bFine = linspace( ...
        min(paramB), ...
        max(paramB), ...
        400);

    yFine = zeros(size(bFine));

    id = compIds{ic};

    funNN = Interpolation.fun{ ...
        id(1),id(2),id(3),id(4)};

    for ib = 1:numel(bFine)

        yFine(ib) = ...
            funNN(bFine(ib),rhoFix);

    end

    plot( ...
        bFine, ...
        yFine, ...
        'r-', ...
        'LineWidth',2);

    xlabel('b');
    ylabel(compNames{ic});

    title(sprintf( ...
        '%s vs b at \\rho = %.4f', ...
        compNames{ic},rhoFix));

    legend( ...
        'HOM', ...
        'NN', ...
        'Location','best');

    grid on
    box on

    exportgraphics( ...
        fig, ...
        fullfile(outFolder, ...
        ['Curve_' compNames{ic} ...
        '_vs_b_rho_050.png']), ...
        'Resolution',300);

    close(fig)

end

fprintf('\nEvaluating NN derivatives dC/db...\n');

fig = figure( ...
    'Color','w', ...
    'Position',[50 50 1400 1000]);

t = tiledlayout(3,2, ...
    'Padding','compact', ...
    'TileSpacing','compact');

for ic = 1:6

    id = compIds{ic};

    dfunNN = Interpolation.dfun{ ...
        id(1),id(2),id(3),id(4)};

    dB_nn = zeros(nB,1);

    for ib = 1:nB

        dNN = dfunNN( ...
            paramB(ib), ...
            rhoFix);

        if iscell(dNN)

            dB_nn(ib) = dNN{1};

        else

            dB_nn(ib) = dNN(1);

        end

    end

    yHom = squeeze(Yhom(irFix,:,ic));

    dB_hom = gradient( ...
        yHom, ...
        paramB);

    nexttile

    plot( ...
        paramB, ...
        dB_hom, ...
        'k-', ...
        'LineWidth',2);

    hold on

    plot( ...
        paramB, ...
        dB_nn, ...
        'r--', ...
        'LineWidth',2);

    xlabel('b');
    ylabel(['d' compNames{ic} '/db']);

    title(compNames{ic});

    legend( ...
        'HOM FD', ...
        'NN', ...
        'Location','best');

    grid on
    box on

end

title(t, ...
    sprintf('dC/db vs b, \\rho = %.4f',rhoFix));

exportgraphics( ...
    fig, ...
    fullfile(outFolder, ...
    'Derivatives_dCdb_vs_b_rho_050.png'), ...
    'Resolution',300);

close(fig)

C1111_hom = squeeze(Yhom(irFix,:,1));
C2222_hom = squeeze(Yhom(irFix,:,4));

ratioHom = C1111_hom./C2222_hom;

fig = figure( ...
    'Color','w', ...
    'Position',[100 100 900 650]);

plot( ...
    paramB, ...
    ratioHom, ...
    'k-', ...
    'LineWidth',2);

hold on

yline(1,'k--','LineWidth',1.5);

xlabel('b');
ylabel('C_{1111}/C_{2222}');

title(sprintf( ...
    'Directional stiffness ratio, \\rho = %.4f', ...
    rhoFix));

grid on
box on

exportgraphics( ...
    fig, ...
    fullfile(outFolder, ...
    'Anisotropy_Ratio_C1111_C2222.png'), ...
    'Resolution',300);

close(fig)

anisotropyHom = ...
    (C1111_hom-C2222_hom) ./ ...
    (0.5*(C1111_hom+C2222_hom));

fig = figure( ...
    'Color','w', ...
    'Position',[100 100 900 650]);

plot( ...
    paramB, ...
    anisotropyHom, ...
    'k-', ...
    'LineWidth',2);

hold on

yline(0,'k--','LineWidth',1.5);

xlabel('b');
ylabel('(C_{1111}-C_{2222}) / [0.5(C_{1111}+C_{2222})]');

title(sprintf( ...
    'Normal anisotropy measure, \\rho = %.4f', ...
    rhoFix));

grid on
box on

exportgraphics( ...
    fig, ...
    fullfile(outFolder, ...
    'Normal_Anisotropy_Measure.png'), ...
    'Resolution',300);

close(fig)

C1212_hom = squeeze(Yhom(irFix,:,6));

[~,ibOne] = min(abs(paramB-1));

C1212ref = C1212_hom(ibOne);

fig = figure( ...
    'Color','w', ...
    'Position',[100 100 900 650]);

plot( ...
    paramB, ...
    C1212_hom/C1212ref, ...
    'k-', ...
    'LineWidth',2);

hold on

yline(1,'k--','LineWidth',1.5);

xlabel('b');
ylabel('C_{1212}/C_{1212}(b \approx 1)');

title(sprintf( ...
    'Relative shear stiffness, \\rho = %.4f', ...
    rhoFix));

grid on
box on

exportgraphics( ...
    fig, ...
    fullfile(outFolder, ...
    'Relative_Shear_Stiffness.png'), ...
    'Resolution',300);

close(fig)

fig = figure( ...
    'Color','w', ...
    'Position',[100 100 900 700]);

imagesc( ...
    paramB, ...
    paramRho, ...
    volFrac);

axis xy

xlabel('b');
ylabel('\rho');

title('Volume fraction');

colorbar

exportgraphics( ...
    fig, ...
    fullfile(outFolder, ...
    'VolumeFraction_Map.png'), ...
    'Resolution',300);

close(fig)

volError = abs(volFrac-paramRho);

fig = figure( ...
    'Color','w', ...
    'Position',[100 100 900 700]);

imagesc( ...
    paramB, ...
    paramRho, ...
    volError);

axis xy

xlabel('b');
ylabel('\rho');

title('|V_f - \rho|');

colorbar

exportgraphics( ...
    fig, ...
    fullfile(outFolder, ...
    'VolumeFraction_Error.png'), ...
    'Resolution',300);

close(fig)

historyFile = 'NNhistory.mat';

if exist(historyFile,'file') == 2

    HH = load(historyFile);

    names = fieldnames(HH);

    fprintf('\nNNhistory.mat variables:\n');

    for i = 1:numel(names)
        fprintf('  %s\n',names{i});
    end

    history = [];

    if isfield(HH,'history')
        history = HH.history;
    end

    if ~isempty(history)

        if isnumeric(history)

            fig = figure( ...
                'Color','w', ...
                'Position',[100 100 900 700]);

            semilogy( ...
                history(:), ...
                'LineWidth',1.5);

            xlabel('Epoch');
            ylabel('Loss');

            title('NN training history');

            grid on
            box on

            exportgraphics( ...
                fig, ...
                fullfile(outFolder, ...
                'NN_Training_History.png'), ...
                'Resolution',300);

            close(fig)

        end

    end

end

fprintf('\n============================================\n');
fprintf('FINISHED\n');
fprintf('============================================\n');
fprintf('All figures saved in:\n');
fprintf('%s\n',outFolder);