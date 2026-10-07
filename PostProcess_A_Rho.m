close all
clc

outFolder = 'PostProcess_A_Rho_Output';

if ~exist(outFolder,'dir')
    mkdir(outFolder);
end

H = load('HomogenizationtwovariablesVA.mat');

Chomog       = H.Chomog;
paramA       = H.paramA(:);
paramRho     = H.paramRho(:);
volFrac      = H.volFrac;
Interpolation = H.Interpolation;

nA   = numel(paramA);
nRho = numel(paramRho);

[Agrid,RHOgrid] = meshgrid(paramA,paramRho);

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

Yhom = zeros(nRho,nA,6);

for ic = 1:6

    id = compIds{ic};

    Yhom(:,:,ic) = squeeze( ...
        Chomog(id(1),id(2),id(3),id(4),:,:) );

end

fprintf('\n============================================\n');
fprintf('POST-PROCESS A-RHO\n');
fprintf('============================================\n');
fprintf('Output folder: %s\n',outFolder);
fprintf('Grid size     : %d rho x %d a\n',nRho,nA);

fprintf('\nEvaluating NN surfaces...\n');

Ynn = zeros(nRho,nA,6);

for ic = 1:6

    id = compIds{ic};

    funNN = Interpolation.fun{ ...
        id(1),id(2),id(3),id(4)};

    for ir = 1:nRho

        for ia = 1:nA

            Ynn(ir,ia,ic) = ...
                funNN(paramA(ia),paramRho(ir));

        end

    end

end

fprintf('NN surfaces evaluated.\n');

for ic = 1:6

    fig = figure( ...
        'Color','w', ...
        'Position',[80 80 900 700]);

    surf( ...
        Agrid, ...
        RHOgrid, ...
        Yhom(:,:,ic), ...
        'EdgeColor','none');

    xlabel('a');
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
        Agrid, ...
        RHOgrid, ...
        Ynn(:,:,ic), ...
        'EdgeColor','none');

    xlabel('a');
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
        Agrid,RHOgrid,Zhom, ...
        'EdgeColor','none');

    xlabel('a');
    ylabel('\rho');
    zlabel(compNames{ic});

    title('HOM');
    view(45,30);
    grid on
    box on
    colorbar

    nexttile

    surf( ...
        Agrid,RHOgrid,Znn, ...
        'EdgeColor','none');

    xlabel('a');
    ylabel('\rho');
    zlabel(compNames{ic});

    title('NN');
    view(45,30);
    grid on
    box on
    colorbar

    nexttile

    imagesc( ...
        paramA, ...
        paramRho, ...
        errAbs);

    axis xy

    xlabel('a');
    ylabel('\rho');

    title('|NN - HOM|');

    colorbar

    nexttile

    imagesc( ...
        paramA, ...
        paramRho, ...
        errGlobal);

    axis xy

    xlabel('a');
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
        paramA, ...
        Yhom(irFix,:,ic), ...
        'k-', ...
        'LineWidth',2);

    hold on

    plot( ...
        paramA, ...
        Ynn(irFix,:,ic), ...
        'r--', ...
        'LineWidth',2);

    xlabel('a');
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
    sprintf('Tensor components vs a, \\rho = %.4f',rhoFix));

exportgraphics( ...
    fig, ...
    fullfile(outFolder, ...
    'Components_vs_a_rho_050.png'), ...
    'Resolution',300);

close(fig)

for ic = 1:6

    fig = figure( ...
        'Color','w', ...
        'Position',[100 100 850 650]);

    plot( ...
        paramA, ...
        Yhom(irFix,:,ic), ...
        'ko-', ...
        'LineWidth',1.5, ...
        'MarkerSize',4);

    hold on

    aFine = linspace( ...
        min(paramA), ...
        max(paramA), ...
        400);

    yFine = zeros(size(aFine));

    id = compIds{ic};

    funNN = Interpolation.fun{ ...
        id(1),id(2),id(3),id(4)};

    for ia = 1:numel(aFine)

        yFine(ia) = ...
            funNN(aFine(ia),rhoFix);

    end

    plot( ...
        aFine, ...
        yFine, ...
        'r-', ...
        'LineWidth',2);

    xlabel('a');
    ylabel(compNames{ic});

    title(sprintf( ...
        '%s vs a at \\rho = %.4f', ...
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
        '_vs_a_rho_050.png']), ...
        'Resolution',300);

    close(fig)

end

fprintf('\nEvaluating NN derivatives dC/da...\n');

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

    dA_nn = zeros(nA,1);

    for ia = 1:nA

        dNN = dfunNN( ...
            paramA(ia), ...
            rhoFix);

        if iscell(dNN)

            dA_nn(ia) = dNN{1};

        else

            dA_nn(ia) = dNN(1);

        end

    end

    yHom = squeeze(Yhom(irFix,:,ic));

    dA_hom = gradient( ...
        yHom, ...
        paramA);

    nexttile

    plot( ...
        paramA, ...
        dA_hom, ...
        'k-', ...
        'LineWidth',2);

    hold on

    plot( ...
        paramA, ...
        dA_nn, ...
        'r--', ...
        'LineWidth',2);

    xlabel('a');
    ylabel(['d' compNames{ic} '/da']);

    title(compNames{ic});

    legend( ...
        'HOM FD', ...
        'NN', ...
        'Location','best');

    grid on
    box on

end

title(t, ...
    sprintf('dC/da vs a, \\rho = %.4f',rhoFix));

exportgraphics( ...
    fig, ...
    fullfile(outFolder, ...
    'Derivatives_dCda_vs_a_rho_050.png'), ...
    'Resolution',300);

close(fig)

fig = figure( ...
    'Color','w', ...
    'Position',[100 100 900 700]);

imagesc( ...
    paramA, ...
    paramRho, ...
    volFrac);

axis xy

xlabel('a');
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
    paramA, ...
    paramRho, ...
    volError);

axis xy

xlabel('a');
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