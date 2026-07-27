%% ================================================================
% POST-PROCESSAMENTO PARA PAPER
%
% Requer no Workspace:
%   ans.f
%   ans.df
%   ans.Chomog
%   ans.paramB
%   ans.paramRho
%   ans.volFrac
%
% Este script NAO limpa o workspace e NAO sobrescreve ans.
% Ele cria uma variavel local H = ans e salva figuras em uma pasta nova.
%% ================================================================

%% ---------------------------------------------------------------
% Configuracoes gerais
%% ---------------------------------------------------------------

H = ans;   % nao sobrescreve ans

timeStamp = datestr(now,'yyyymmdd_HHMMSS');
outDir = ['PaperPlots_', timeStamp];

if ~exist(outDir,'dir')
    mkdir(outDir);
end

saveFigures = true;
figFormat   = 'png';   % pode trocar para 'pdf' ou 'eps'

fprintf('\n============================================================\n');
fprintf('POST-PROCESSAMENTO PARA PAPER\n');
fprintf('============================================================\n');
fprintf('Pasta de saida: %s\n', outDir);

paramB   = H.paramB(:)';
paramRho = H.paramRho(:)';

bMin   = min(paramB);
bMax   = max(paramB);
rhoMin = min(paramRho);
rhoMax = max(paramRho);

[Bgrid,Rgrid] = meshgrid(paramB,paramRho);

compMap = {
    [1,1,1,1], 'C1111';
    [2,2,2,2], 'C2222';
    [1,1,2,2], 'C1122';
    [1,2,1,2], 'C1212';
    [1,1,1,2], 'C1112';
    [2,2,1,2], 'C2212'
};

nComp = size(compMap,1);
nRho  = numel(paramRho);
nB    = numel(paramB);

%% ---------------------------------------------------------------
% Extrair C FEM e C NN
%% ---------------------------------------------------------------

Cref  = zeros(nRho,nB,nComp);
Cpred = zeros(nRho,nB,nComp);

for m = 1:nComp
    idx  = compMap{m,1};
    fH   = H.f{idx(1),idx(2),idx(3),idx(4)};

    Cref(:,:,m) = squeeze(H.Chomog(idx(1),idx(2),idx(3),idx(4),:,:));
    Cpred(:,:,m) = fH(Bgrid,Rgrid);

    if ~isequal(size(Cref(:,:,m)),size(Cpred(:,:,m)))
        error('Dimensoes inconsistentes para %s.',compMap{m,2});
    end
end

CerrAbs = abs(Cpred-Cref);

CerrRel = zeros(size(CerrAbs));
for m = 1:nComp
    scale = max(abs(Cref(:,:,m)),[],'all');
    CerrRel(:,:,m) = CerrAbs(:,:,m)./max(abs(Cref(:,:,m)),1e-8*scale + 1e-12);
end

%% ---------------------------------------------------------------
% Tabela de erros das componentes C
%% ---------------------------------------------------------------

fprintf('\n============================================================\n');
fprintf('ERROS DAS COMPONENTES DE C\n');
fprintf('============================================================\n');

errorTable = table();

for m = 1:nComp
    name = compMap{m,2};

    e = Cpred(:,:,m)-Cref(:,:,m);
    rmse = sqrt(mean(e(:).^2));
    mae  = mean(abs(e(:)));
    maxAbs = max(abs(e(:)));

    refRMS = sqrt(mean(Cref(:,:,m).^2,'all'));
    nrmseRMS = rmse/max(refRMS,1e-12);

    refRange = max(Cref(:,:,m),[],'all') - min(Cref(:,:,m),[],'all');
    nrmseRange = rmse/max(refRange,1e-12);

    r2Den = sum((Cref(:,:,m)-mean(Cref(:,:,m),'all')).^2,'all');
    r2Num = sum(e.^2,'all');

    if r2Den > 0
        R2 = 1 - r2Num/r2Den;
    else
        R2 = NaN;
    end

    errorTable = [errorTable; table( ...
        string(name), rmse, mae, maxAbs, nrmseRMS, nrmseRange, R2, ...
        'VariableNames',{'Component','RMSE','MAE','MaxAbs','NRMSE_RMS','NRMSE_Range','R2'})]; %#ok<AGROW>

    fprintf('%s: RMSE = %.4e, NRMSE-RMS = %.4e, MaxAbs = %.4e, R2 = %.8f\n', ...
        name, rmse, nrmseRMS, maxAbs, R2);
end

disp(errorTable);

%% ---------------------------------------------------------------
% Plotar superficies C FEM, C NN e erro
%% ---------------------------------------------------------------

for m = 1:nComp
    name = compMap{m,2};

    plotSurfaceMap(Bgrid,Rgrid,Cref(:,:,m), ...
        ['$' name '^{FEM}$'], '$b$', '$\rho$');
    saveCurrentFigure(saveFigures,outDir,['Surface_FEM_',name],figFormat);

    plotSurfaceMap(Bgrid,Rgrid,Cpred(:,:,m), ...
        ['$' name '^{NN}$'], '$b$', '$\rho$');
    saveCurrentFigure(saveFigures,outDir,['Surface_NN_',name],figFormat);

    plotSurfaceMap(Bgrid,Rgrid,CerrAbs(:,:,m), ...
        ['Absolute error ', name], '$b$', '$\rho$');
    saveCurrentFigure(saveFigures,outDir,['ErrorAbs_',name],figFormat);

    plotSurfaceMap(Bgrid,Rgrid,CerrRel(:,:,m), ...
        ['Relative error ', name], '$b$', '$\rho$');
    saveCurrentFigure(saveFigures,outDir,['ErrorRel_',name],figFormat);
end

%% ---------------------------------------------------------------
% Plot combinado das 6 componentes C FEM e NN
%% ---------------------------------------------------------------

figure('Name','All C components FEM','Color','w');
tiledlayout(2,3,'Padding','compact','TileSpacing','compact');

for m = 1:nComp
    nexttile;
    surf(Bgrid,Rgrid,Cref(:,:,m),'EdgeColor','none');
    view(35,25);
    xlabel('$b$','Interpreter','latex');
    ylabel('$\rho$','Interpreter','latex');
    zlabel(compMap{m,2},'Interpreter','none');
    title([compMap{m,2}, ' FEM'],'Interpreter','none');
    colorbar;
end
saveCurrentFigure(saveFigures,outDir,'All_C_FEM',figFormat);

figure('Name','All C components NN','Color','w');
tiledlayout(2,3,'Padding','compact','TileSpacing','compact');

for m = 1:nComp
    nexttile;
    surf(Bgrid,Rgrid,Cpred(:,:,m),'EdgeColor','none');
    view(35,25);
    xlabel('$b$','Interpreter','latex');
    ylabel('$\rho$','Interpreter','latex');
    zlabel(compMap{m,2},'Interpreter','none');
    title([compMap{m,2}, ' NN'],'Interpreter','none');
    colorbar;
end
saveCurrentFigure(saveFigures,outDir,'All_C_NN',figFormat);

%% ---------------------------------------------------------------
% Reconstruir Cholesky L para FEM e NN
%% ---------------------------------------------------------------

Lnames = {'L11','L21','L22','L31','L32','L33'};

Lref  = zeros(nRho,nB,6);
Lpred = zeros(nRho,nB,6);

minEigRef  = Inf;
minEigPred = Inf;
nFailRef   = 0;
nFailPred  = 0;

for ir = 1:nRho
    for ib = 1:nB

        CmatRef = voigtMatrixFromComponents(squeeze(Cref(ir,ib,:)));
        CmatPred = voigtMatrixFromComponents(squeeze(Cpred(ir,ib,:)));

        CmatRef = 0.5*(CmatRef+CmatRef');
        CmatPred = 0.5*(CmatPred+CmatPred');

        evRef = eig(CmatRef);
        evPred = eig(CmatPred);

        minEigRef  = min(minEigRef,min(evRef));
        minEigPred = min(minEigPred,min(evPred));

        [Lr,flagR] = cholSPD(CmatRef);
        [Lp,flagP] = cholSPD(CmatPred);

        if flagR ~= 0
            nFailRef = nFailRef+1;
        end

        if flagP ~= 0
            nFailPred = nFailPred+1;
        end

        Lref(ir,ib,:) = [Lr(1,1), Lr(2,1), Lr(2,2), Lr(3,1), Lr(3,2), Lr(3,3)];
        Lpred(ir,ib,:) = [Lp(1,1), Lp(2,1), Lp(2,2), Lp(3,1), Lp(3,2), Lp(3,3)];
    end
end

fprintf('\n============================================================\n');
fprintf('POSITIVIDADE DEFINIDA / CHOLESKY\n');
fprintf('============================================================\n');
fprintf('Menor autovalor FEM: %.6e\n',minEigRef);
fprintf('Menor autovalor NN:  %.6e\n',minEigPred);
fprintf('Falhas Cholesky FEM: %d\n',nFailRef);
fprintf('Falhas Cholesky NN:  %d\n',nFailPred);

%% ---------------------------------------------------------------
% Plotar componentes de L
%% ---------------------------------------------------------------

for q = 1:6
    plotSurfaceMap(Bgrid,Rgrid,Lref(:,:,q), ...
        ['$' Lnames{q} '^{FEM}$'], '$b$', '$\rho$');
    saveCurrentFigure(saveFigures,outDir,['L_FEM_',Lnames{q}],figFormat);

    plotSurfaceMap(Bgrid,Rgrid,Lpred(:,:,q), ...
        ['$' Lnames{q} '^{NN}$'], '$b$', '$\rho$');
    saveCurrentFigure(saveFigures,outDir,['L_NN_',Lnames{q}],figFormat);

    plotSurfaceMap(Bgrid,Rgrid,abs(Lpred(:,:,q)-Lref(:,:,q)), ...
        ['Absolute error ', Lnames{q}], '$b$', '$\rho$');
    saveCurrentFigure(saveFigures,outDir,['L_ErrorAbs_',Lnames{q}],figFormat);
end

figure('Name','All L components NN','Color','w');
tiledlayout(2,3,'Padding','compact','TileSpacing','compact');

for q = 1:6
    nexttile;
    surf(Bgrid,Rgrid,Lpred(:,:,q),'EdgeColor','none');
    view(35,25);
    xlabel('$b$','Interpreter','latex');
    ylabel('$\rho$','Interpreter','latex');
    zlabel(Lnames{q},'Interpreter','none');
    title([Lnames{q}, ' NN'],'Interpreter','none');
    colorbar;
end
saveCurrentFigure(saveFigures,outDir,'All_L_NN',figFormat);

%% ---------------------------------------------------------------
% Simetria do tensor Chomog original
%% ---------------------------------------------------------------

fprintf('\n============================================================\n');
fprintf('SIMETRIA DO TENSOR HOMOGENEIZADO ORIGINAL\n');
fprintf('============================================================\n');

maxMinorIJ = 0;
maxMinorKL = 0;
maxMajor   = 0;
maxNormC   = 0;

for ir = 1:nRho
    for ib = 1:nB
        Ct = H.Chomog(:,:,:,:,ir,ib);
        maxNormC = max(maxNormC,norm(Ct(:)));

        errIJ = 0;
        errKL = 0;
        errMaj = 0;

        for i = 1:2
            for j = 1:2
                for k = 1:2
                    for l = 1:2
                        errIJ  = max(errIJ, abs(Ct(i,j,k,l)-Ct(j,i,k,l)));
                        errKL  = max(errKL, abs(Ct(i,j,k,l)-Ct(i,j,l,k)));
                        errMaj = max(errMaj,abs(Ct(i,j,k,l)-Ct(k,l,i,j)));
                    end
                end
            end
        end

        maxMinorIJ = max(maxMinorIJ,errIJ);
        maxMinorKL = max(maxMinorKL,errKL);
        maxMajor   = max(maxMajor,errMaj);
    end
end

fprintf('Max minor symmetry error ij: %.6e\n',maxMinorIJ);
fprintf('Max minor symmetry error kl: %.6e\n',maxMinorKL);
fprintf('Max major symmetry error:    %.6e\n',maxMajor);
fprintf('Relative major symmetry:     %.6e\n',maxMajor/max(maxNormC,1e-12));

%% ---------------------------------------------------------------
% Anisotropia, Zener-like ratio e indicadores
%% ---------------------------------------------------------------

C11 = Cpred(:,:,1);
C22 = Cpred(:,:,2);
C12 = Cpred(:,:,3);
C66 = Cpred(:,:,4);
C16 = Cpred(:,:,5);
C26 = Cpred(:,:,6);

anisRatio = C11./max(C22,1e-14);

% Zener-like ratio para 2D:
% Para isotropia: C11 = C22, C66 = (C11-C12)/2, logo Z = 1.
zenerLike = 2*C66./max(C11-C12,1e-14);

% Indice simples de anisotropia:
% mede desvio em relacao a relacoes isotropicas principais.
isoDev = sqrt( ...
    (C11-C22).^2 + ...
    (2*C66 - (C11-C12)).^2 + ...
    C16.^2 + C26.^2 );

isoScale = sqrt(C11.^2 + C22.^2 + C12.^2 + C66.^2 + C16.^2 + C26.^2);
anisIndex = isoDev./max(isoScale,1e-14);

plotSurfaceMap(Bgrid,Rgrid,anisRatio, ...
    '$C_{1111}/C_{2222}$', '$b$', '$\rho$');
saveCurrentFigure(saveFigures,outDir,'Anisotropy_C11_over_C22',figFormat);

plotSurfaceMap(Bgrid,Rgrid,zenerLike, ...
    '$Z_{2D}=2C_{1212}/(C_{1111}-C_{1122})$', '$b$', '$\rho$');
saveCurrentFigure(saveFigures,outDir,'ZenerLike_2D',figFormat);

plotSurfaceMap(Bgrid,Rgrid,anisIndex, ...
    'Normalized anisotropy index', '$b$', '$\rho$');
saveCurrentFigure(saveFigures,outDir,'AnisotropyIndex',figFormat);

fprintf('\n============================================================\n');
fprintf('ANISOTROPIA\n');
fprintf('============================================================\n');
fprintf('C11/C22: min %.4e, max %.4e\n',min(anisRatio(:)),max(anisRatio(:)));
fprintf('Zener-like: min %.4e, max %.4e\n',min(zenerLike(:)),max(zenerLike(:)));
fprintf('Anisotropy index: min %.4e, max %.4e\n',min(anisIndex(:)),max(anisIndex(:)));

%% ---------------------------------------------------------------
% Curvas para rho fixo: C, Zener, anisotropia vs b
%% ---------------------------------------------------------------

rhoCuts = [0.2, 0.4, 0.6, 0.8];

for rr = 1:numel(rhoCuts)

    [~,ir] = min(abs(paramRho-rhoCuts(rr)));
    rhoVal = paramRho(ir);

    figure('Name',['C_vs_b_rho_',num2str(rhoVal)],'Color','w');
    hold on; grid on; box on;

    plot(paramB,squeeze(Cpred(ir,:,1)),'LineWidth',1.8);
    plot(paramB,squeeze(Cpred(ir,:,2)),'LineWidth',1.8);
    plot(paramB,squeeze(Cpred(ir,:,3)),'LineWidth',1.8);
    plot(paramB,squeeze(Cpred(ir,:,4)),'LineWidth',1.8);
    plot(paramB,squeeze(Cpred(ir,:,5)),'LineWidth',1.8);
    plot(paramB,squeeze(Cpred(ir,:,6)),'LineWidth',1.8);

    xlabel('$b$','Interpreter','latex');
    ylabel('$C^H$ component','Interpreter','latex');
    title(sprintf('Homogenized tensor components at rho = %.3f',rhoVal),'Interpreter','none');
    legend({'C1111','C2222','C1122','C1212','C1112','C2212'}, ...
        'Location','best','Interpreter','none');

    saveCurrentFigure(saveFigures,outDir,sprintf('C_vs_b_rho_%0.3f',rhoVal),figFormat);

    figure('Name',['Anisotropy_vs_b_rho_',num2str(rhoVal)],'Color','w');
    hold on; grid on; box on;

    plot(paramB,anisRatio(ir,:),'LineWidth',2);
    plot(paramB,zenerLike(ir,:),'LineWidth',2);
    plot(paramB,anisIndex(ir,:),'LineWidth',2);

    xlabel('$b$','Interpreter','latex');
    ylabel('Indicator','Interpreter','latex');
    title(sprintf('Anisotropy indicators at rho = %.3f',rhoVal),'Interpreter','none');
    legend({'C1111/C2222','Zener-like','Anisotropy index'}, ...
        'Location','best','Interpreter','none');

    saveCurrentFigure(saveFigures,outDir,sprintf('Anisotropy_vs_b_rho_%0.3f',rhoVal),figFormat);
end

%% ---------------------------------------------------------------
% Curvas para b fixo: C, Zener, anisotropia vs rho
%% ---------------------------------------------------------------

bCuts = [-0.6, -0.3, 0, 0.3, 0.6];

for bb = 1:numel(bCuts)

    [~,ib] = min(abs(paramB-bCuts(bb)));
    bVal = paramB(ib);

    figure('Name',['C_vs_rho_b_',num2str(bVal)],'Color','w');
    hold on; grid on; box on;

    plot(paramRho,squeeze(Cpred(:,ib,1)),'LineWidth',1.8);
    plot(paramRho,squeeze(Cpred(:,ib,2)),'LineWidth',1.8);
    plot(paramRho,squeeze(Cpred(:,ib,3)),'LineWidth',1.8);
    plot(paramRho,squeeze(Cpred(:,ib,4)),'LineWidth',1.8);

    xlabel('$\rho$','Interpreter','latex');
    ylabel('$C^H$ component','Interpreter','latex');
    title(sprintf('Main tensor components at b = %.3f',bVal),'Interpreter','none');
    legend({'C1111','C2222','C1122','C1212'}, ...
        'Location','best','Interpreter','none');

    saveCurrentFigure(saveFigures,outDir,sprintf('C_vs_rho_b_%0.3f',bVal),figFormat);

    figure('Name',['Anisotropy_vs_rho_b_',num2str(bVal)],'Color','w');
    hold on; grid on; box on;

    plot(paramRho,anisRatio(:,ib),'LineWidth',2);
    plot(paramRho,zenerLike(:,ib),'LineWidth',2);
    plot(paramRho,anisIndex(:,ib),'LineWidth',2);

    xlabel('$\rho$','Interpreter','latex');
    ylabel('Indicator','Interpreter','latex');
    title(sprintf('Anisotropy indicators at b = %.3f',bVal),'Interpreter','none');
    legend({'C1111/C2222','Zener-like','Anisotropy index'}, ...
        'Location','best','Interpreter','none');

    saveCurrentFigure(saveFigures,outDir,sprintf('Anisotropy_vs_rho_b_%0.3f',bVal),figFormat);
end

%% ---------------------------------------------------------------
% Periodicidade geometrica: tiling de celulas transformadas
%% ---------------------------------------------------------------

rhoGeom = [0.3, 0.6, 0.85];
bGeom   = [-0.6, 0, 0.6];

for ibg = 1:numel(bGeom)
    for irg = 1:numel(rhoGeom)

        figure('Name','Periodic cell tiling','Color','w');
        hold on; axis equal; box on; grid on;

        plotPeriodicTiling(bGeom(ibg),rhoGeom(irg),4,4);

        title(sprintf('Periodic tiling: b = %.2f, rho = %.2f', ...
            bGeom(ibg),rhoGeom(irg)), ...
            'Interpreter','none');

        xlabel('$x_1$','Interpreter','latex');
        ylabel('$x_2$','Interpreter','latex');

        saveCurrentFigure(saveFigures,outDir, ...
            sprintf('PeriodicTiling_b_%0.2f_rho_%0.2f',bGeom(ibg),rhoGeom(irg)),figFormat);
    end
end

%% ---------------------------------------------------------------
% Salvar tabelas de metricas em arquivos separados
%% ---------------------------------------------------------------

tableFile = fullfile(outDir,'ErrorTable_C_components.csv');
writetable(errorTable,tableFile);

metricFile = fullfile(outDir,'AnisotropyMetrics.mat');
postResults = struct();
postResults.errorTable = errorTable;
postResults.minEigRef = minEigRef;
postResults.minEigPred = minEigPred;
postResults.maxMinorIJ = maxMinorIJ;
postResults.maxMinorKL = maxMinorKL;
postResults.maxMajor = maxMajor;
postResults.anisRatio = anisRatio;
postResults.zenerLike = zenerLike;
postResults.anisIndex = anisIndex;
postResults.paramB = paramB;
postResults.paramRho = paramRho;

save(metricFile,'postResults');

fprintf('\n============================================================\n');
fprintf('POST-PROCESSAMENTO FINALIZADO\n');
fprintf('Figuras e tabelas em: %s\n', outDir);
fprintf('Nao foi feito clear e ans nao foi sobrescrito.\n');
fprintf('============================================================\n');

%% ================================================================
% FUNCOES LOCAIS
%% ================================================================

function Cmat = voigtMatrixFromComponents(c)
    C11 = c(1);
    C22 = c(2);
    C12 = c(3);
    C66 = c(4);
    C16 = c(5);
    C26 = c(6);

    Cmat = [
        C11, C12, C16;
        C12, C22, C26;
        C16, C26, C66
    ];
end

function [L,flag] = cholSPD(Cmat)
    Cmat = 0.5*(Cmat+Cmat');
    [L,flag] = chol(Cmat,'lower');

    if flag ~= 0
        jitter = max(1e-12,1e-10*norm(Cmat,'fro'));
        [L,flag] = chol(Cmat + jitter*eye(size(Cmat)),'lower');
    end

    if flag ~= 0
        [V,D] = eig(Cmat);
        d = max(diag(D),1e-14);
        Cfix = V*diag(d)*V';
        Cfix = 0.5*(Cfix+Cfix');
        [L,flag] = chol(Cfix,'lower');
    end
end

function plotSurfaceMap(X,Y,Z,plotTitle,xlab,ylab)
    figure('Color','w');
    surf(X,Y,Z,'EdgeColor','none');
    view(35,25);
    xlabel(xlab,'Interpreter','latex');
    ylabel(ylab,'Interpreter','latex');
    zlabel(plotTitle,'Interpreter','latex');
    title(plotTitle,'Interpreter','latex');
    colorbar;
    grid on;
    box on;
end

function saveCurrentFigure(saveFigures,outDir,name,figFormat)
    if saveFigures
        fileName = fullfile(outDir,[name,'.',figFormat]);
        exportgraphics(gcf,fileName,'Resolution',300);
    end
end

function plotPeriodicTiling(b,rho,n1,n2)

    a = exp(b.^2);
    d = (1+b.^2)./a;

    v1 = [a; b];
    v2 = [b; d];

    rhoClipped = min(max(real(rho),1e-6),1-1e-6);
    holeHalfSize = sqrt(1-rhoClipped)/2;

    % quadrado da celula em coordenadas locais xi
    cellXi = [0 1 1 0 0; 0 0 1 1 0];

    % furo quadrado em coordenadas locais centradas
    h = holeHalfSize;
    holeXi = [
        0.5-h, 0.5+h, 0.5+h, 0.5-h, 0.5-h;
        0.5-h, 0.5-h, 0.5+h, 0.5+h, 0.5-h
    ];

    T = [v1, v2];

    for i = 0:n1-1
        for j = 0:n2-1

            shift = i*v1 + j*v2;

            cellXY = T*cellXi + shift;
            holeXY = T*holeXi + shift;

            plot(cellXY(1,:),cellXY(2,:),'k-','LineWidth',0.8);
            plot(holeXY(1,:),holeXY(2,:),'r-','LineWidth',1.2);

        end
    end

    quiver(0,0,v1(1),v1(2),0,'LineWidth',2,'MaxHeadSize',0.5);
    quiver(0,0,v2(1),v2(2),0,'LineWidth',2,'MaxHeadSize',0.5);

    text(v1(1),v1(2),'  $a_1(b)$','Interpreter','latex','FontSize',12);
    text(v2(1),v2(2),'  $a_2(b)$','Interpreter','latex','FontSize',12);

    axis equal;
end