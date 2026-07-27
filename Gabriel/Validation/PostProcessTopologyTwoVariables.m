%% ================================================================
% POST-PROCESSAMENTO DA OTIMIZACAO TOPOLOGICA b-rho
%
% Requer no workspace:
%   TH = resultado da otimizacao
%
% Este script NAO limpa o workspace e NAO sobrescreve ans.
%% ================================================================

%% ---------------------------------------------------------------
% Segurança
%% ---------------------------------------------------------------

if ~exist('TH','var')
    error('Crie primeiro TH = ans; antes de rodar este script.');
end

timeStamp = datestr(now,'yyyymmdd_HHMMSS');
outDir = ['TopologyPlots_',timeStamp];

if ~exist(outDir,'dir')
    mkdir(outDir);
end

fprintf('\n============================================================\n');
fprintf('POST-PROCESSAMENTO DA TOPOLOGIA\n');
fprintf('============================================================\n');
fprintf('Pasta de saida: %s\n',outDir);

%% ---------------------------------------------------------------
% Extrair campos b e rho
%% ---------------------------------------------------------------

bFun   = TH.designVariable.funB;
rhoFun = TH.designVariable.funRho;

bVals   = bFun.fValues;
rhoVals = rhoFun.fValues;

fprintf('\n============================================================\n');
fprintf('RANGES DOS CAMPOS OTIMIZADOS\n');
fprintf('============================================================\n');

fprintf('b   min/max  = %.6f %.6f\n', min(bVals), max(bVals));
fprintf('rho min/max  = %.6f %.6f\n', min(rhoVals), max(rhoVals));
fprintf('mean rho     = %.6f\n', mean(rhoVals));
fprintf('std rho      = %.6f\n', std(rhoVals));
fprintf('mean b       = %.6f\n', mean(bVals));
fprintf('std b        = %.6f\n', std(bVals));

fprintf('\nPercentis:\n');
fprintf('b   p01/p50/p99  = %.6f %.6f %.6f\n', ...
    prctile(bVals,1), prctile(bVals,50), prctile(bVals,99));
fprintf('rho p01/p50/p99  = %.6f %.6f %.6f\n', ...
    prctile(rhoVals,1), prctile(rhoVals,50), prctile(rhoVals,99));

%% ---------------------------------------------------------------
% Plot usando o metodo interno .plot()
%% ---------------------------------------------------------------

figure('Color','w');
bFun.plot();
title('Optimized b field','Interpreter','latex');
colorbar;
exportgraphics(gcf,fullfile(outDir,'b_field_internal_plot.png'),'Resolution',300);

figure('Color','w');
rhoFun.plot();
title('Optimized $\rho$ field','Interpreter','latex');
colorbar;
exportgraphics(gcf,fullfile(outDir,'rho_field_internal_plot.png'),'Resolution',300);

%% ---------------------------------------------------------------
% Plots coloridos e preto/branco usando malha
%% ---------------------------------------------------------------

mesh = bFun.mesh;

if isprop(mesh,'coord') || isfield(mesh,'coord')
    coord = mesh.coord;
else
    coord = [];
end

% Tentativa 1: usar coordenadas dos nós se baterem com fValues
if ~isempty(coord) && size(coord,1) == numel(bVals)

    x = coord(:,1);
    y = coord(:,2);

    figure('Color','w');
    scatter(x,y,20,bVals,'filled');
    axis equal tight;
    title('Optimized b field','Interpreter','latex');
    xlabel('$x$','Interpreter','latex');
    ylabel('$y$','Interpreter','latex');
    colorbar;
    exportgraphics(gcf,fullfile(outDir,'b_field_color_scatter.png'),'Resolution',300);

    figure('Color','w');
    scatter(x,y,20,rhoVals,'filled');
    axis equal tight;
    title('Optimized $\rho$ field','Interpreter','latex');
    xlabel('$x$','Interpreter','latex');
    ylabel('$y$','Interpreter','latex');
    colorbar;
    exportgraphics(gcf,fullfile(outDir,'rho_field_color_scatter.png'),'Resolution',300);

    figure('Color','w');
    scatter(x,y,20,bVals,'filled');
    axis equal tight;
    colormap(gray);
    title('Optimized b field -- grayscale','Interpreter','latex');
    xlabel('$x$','Interpreter','latex');
    ylabel('$y$','Interpreter','latex');
    colorbar;
    exportgraphics(gcf,fullfile(outDir,'b_field_grayscale_scatter.png'),'Resolution',300);

    figure('Color','w');
    scatter(x,y,20,rhoVals,'filled');
    axis equal tight;
    colormap(gray);
    title('Optimized $\rho$ field -- grayscale','Interpreter','latex');
    xlabel('$x$','Interpreter','latex');
    ylabel('$y$','Interpreter','latex');
    colorbar;
    exportgraphics(gcf,fullfile(outDir,'rho_field_grayscale_scatter.png'),'Resolution',300);

else
    fprintf('\nNao foi possivel fazer scatter por coordenadas nodais.\n');
    fprintf('Os plots internos .plot() foram gerados corretamente.\n');
end

%% ---------------------------------------------------------------
% Histograma dos campos
%% ---------------------------------------------------------------

figure('Color','w');
histogram(bVals,40);
grid on; box on;
xlabel('$b$','Interpreter','latex');
ylabel('count','Interpreter','latex');
title('Histogram of optimized b','Interpreter','latex');
exportgraphics(gcf,fullfile(outDir,'histogram_b.png'),'Resolution',300);

figure('Color','w');
histogram(rhoVals,40);
grid on; box on;
xlabel('$\rho$','Interpreter','latex');
ylabel('count','Interpreter','latex');
title('Histogram of optimized $\rho$','Interpreter','latex');
exportgraphics(gcf,fullfile(outDir,'histogram_rho.png'),'Resolution',300);

%% ---------------------------------------------------------------
% Checagem de bounds
%% ---------------------------------------------------------------

bLower   = -0.8;
bUpper   =  0.8;
rhoLower =  0.001;
rhoUpper =  0.95;

tol = 1e-8;

nBelowB   = sum(bVals   < bLower   - tol);
nAboveB   = sum(bVals   > bUpper   + tol);
nBelowRho = sum(rhoVals < rhoLower - tol);
nAboveRho = sum(rhoVals > rhoUpper + tol);

fprintf('\n============================================================\n');
fprintf('CHECAGEM DE BOUNDS\n');
fprintf('============================================================\n');

fprintf('b abaixo do limite:    %d\n', nBelowB);
fprintf('b acima do limite:     %d\n', nAboveB);
fprintf('rho abaixo do limite:  %d\n', nBelowRho);
fprintf('rho acima do limite:   %d\n', nAboveRho);

if nBelowB+nAboveB+nBelowRho+nAboveRho == 0
    fprintf('Bounds OK.\n');
else
    fprintf('Atencao: algum campo passou dos bounds.\n');
end

%% ---------------------------------------------------------------
% Volume aproximado
%% ---------------------------------------------------------------

fprintf('\n============================================================\n');
fprintf('VOLUME / DENSIDADE MEDIA\n');
fprintf('============================================================\n');

fprintf('mean(rho) = %.6f\n', mean(rhoVals));

% Se a malha tiver volumes de elemento, tente calcular media ponderada
if isprop(mesh,'volume') || isfield(mesh,'volume')
    try
        w = mesh.volume(:);
        if numel(w) == numel(rhoVals)
            volWeighted = sum(rhoVals(:).*w)/sum(w);
            fprintf('volume ponderado por mesh.volume = %.6f\n',volWeighted);
        end
    catch
    end
end

%% ---------------------------------------------------------------
% Salvar resumo num .mat e .txt sem sobrescrever ans
%% ---------------------------------------------------------------

topologyPost = struct();

topologyPost.bMin   = min(bVals);
topologyPost.bMax   = max(bVals);
topologyPost.rhoMin = min(rhoVals);
topologyPost.rhoMax = max(rhoVals);
topologyPost.meanRho = mean(rhoVals);
topologyPost.stdRho  = std(rhoVals);
topologyPost.meanB   = mean(bVals);
topologyPost.stdB    = std(bVals);
topologyPost.nBelowB = nBelowB;
topologyPost.nAboveB = nAboveB;
topologyPost.nBelowRho = nBelowRho;
topologyPost.nAboveRho = nAboveRho;

save(fullfile(outDir,'TopologyPostResults.mat'),'topologyPost');

txtFile = fullfile(outDir,'TopologySummary.txt');
fid = fopen(txtFile,'w');

fprintf(fid,'POST-PROCESSAMENTO DA TOPOLOGIA\n');
fprintf(fid,'Data: %s\n',datestr(now));
fprintf(fid,'b min/max   = %.10f %.10f\n',topologyPost.bMin,topologyPost.bMax);
fprintf(fid,'rho min/max = %.10f %.10f\n',topologyPost.rhoMin,topologyPost.rhoMax);
fprintf(fid,'mean rho    = %.10f\n',topologyPost.meanRho);
fprintf(fid,'std rho     = %.10f\n',topologyPost.stdRho);
fprintf(fid,'mean b      = %.10f\n',topologyPost.meanB);
fprintf(fid,'std b       = %.10f\n',topologyPost.stdB);
fprintf(fid,'nBelowB     = %d\n',topologyPost.nBelowB);
fprintf(fid,'nAboveB     = %d\n',topologyPost.nAboveB);
fprintf(fid,'nBelowRho   = %d\n',topologyPost.nBelowRho);
fprintf(fid,'nAboveRho   = %d\n',topologyPost.nAboveRho);

fclose(fid);

fprintf('\n============================================================\n');
fprintf('POST-PROCESSAMENTO FINALIZADO\n');
fprintf('Resultados em: %s\n',outDir);
fprintf('ans nao foi sobrescrito.\n');
fprintf('============================================================\n');