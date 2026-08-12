%% RunDehomogenizationSavedCases
% Reconstrucao/dehomogeneizacao dos tres casos:
%   1) B-only        : b variavel, rho = 0.4
%   2) Rho-only      : b = 0, rho variavel
%   3) Two variables : b e rho variaveis
%
% O script usa somente os arquivos .mat salvos em
% TopologyComparisonResults. Nao executa novamente a otimizacao.
%
% Saidas para cada caso:
%   - arquivo .mat com os campos reconstruidos;
%   - figura .png e .fig da geometria solida;
%   - figura .png e .fig do level set;
%   - arquivo .vti para abrir no ParaView.
%
% No ParaView:
%   1) abra o arquivo .vti;
%   2) selecione o campo "Solid" para visualizar a geometria;
%   3) use "LevelSet" e o filtro Contour com valor 0 para a fronteira.

clear;
close all;
clc;

%% ============================================================
%  CONFIGURACAO
%  ============================================================

resultsFolder = fullfile(pwd,'TopologyComparisonResults');
outputFolder  = fullfile(pwd,'DehomogenizationResults');

if ~exist(resultsFolder,'dir')
    error('Pasta nao encontrada: %s',resultsFolder);
end

if ~exist(outputFolder,'dir')
    mkdir(outputFolder);
end

% Tamanho da celula microscopica no dominio macro.
% Comece com 0.10. Depois pode reduzir para 0.07.
epsMicro = 0.10;

% Grade regular usada na reconstrucao e no VTI.
% Valores moderados para reduzir memoria e tempo.
Nx = 500;
Ny = 250;

% Suavizacao leve do campo b antes da reconstrucao.
nSmoothB = 1;

fprintf('\n');
fprintf('============================================================\n');
fprintf('DEHOMOGENIZATION OF SAVED OPTIMIZED FIELDS\n');
fprintf('============================================================\n');
fprintf('epsMicro = %.6f\n',epsMicro);
fprintf('grid     = %d x %d\n',Nx,Ny);
fprintf('output   = %s\n',outputFolder);
fprintf('============================================================\n');

%% ============================================================
%  CARREGAR RESULTADOS
%  ============================================================

load(fullfile(resultsFolder,'Result_BOnly.mat'), ...
    'resultBOnly');

load(fullfile(resultsFolder,'Result_RhoOnly.mat'), ...
    'resultRhoOnly');

load(fullfile(resultsFolder,'Result_TwoVariable.mat'), ...
    'resultTwoVariable');

%% ============================================================
%  CASO 1: B-ONLY
%  ============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('CASE 1/3 - B-ONLY\n');
fprintf('============================================================\n');

meshB = resultBOnly.mesh;

% Usa o parametro geometrico real reconstruido a partir de beta.
bB = resultBOnly.bGeomSmooth(:);

% Densidade fixa.
rhoB = 0.4*ones(size(bB));

dataB = reconstructDehomogenizedField( ...
    meshB,bB,rhoB,epsMicro,Nx,Ny,nSmoothB,false);

saveCaseResults( ...
    dataB,outputFolder,'Dehomog_BOnly',epsMicro);

clear dataB meshB bB rhoB

%% ============================================================
%  CASO 2: RHO-ONLY
%  ============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('CASE 2/3 - RHO-ONLY\n');
fprintf('============================================================\n');

meshRho = resultRhoOnly.mesh;

rhoRho = resultRhoOnly.rhoSmooth(:);
bRho   = zeros(size(rhoRho));

% Como b = 0 em todo o dominio, as fases sao exatas e nao
% precisam de LSQR.
dataRho = reconstructDehomogenizedField( ...
    meshRho,bRho,rhoRho,epsMicro,Nx,Ny,0,true);

saveCaseResults( ...
    dataRho,outputFolder,'Dehomog_RhoOnly',epsMicro);

clear dataRho meshRho bRho rhoRho

%% ============================================================
%  CASO 3: TWO VARIABLES
%  ============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('CASE 3/3 - TWO VARIABLES\n');
fprintf('============================================================\n');

meshTwo = resultTwoVariable.mesh;

bTwo   = resultTwoVariable.bSmooth(:);
rhoTwo = resultTwoVariable.rhoSmooth(:);

dataTwo = reconstructDehomogenizedField( ...
    meshTwo,bTwo,rhoTwo,epsMicro,Nx,Ny,nSmoothB,false);

saveCaseResults( ...
    dataTwo,outputFolder,'Dehomog_TwoVariable',epsMicro);

clear dataTwo meshTwo bTwo rhoTwo

fprintf('\n');
fprintf('============================================================\n');
fprintf('ALL DEHOMOGENIZATION CASES COMPLETED\n');
fprintf('Results saved in:\n%s\n',outputFolder);
fprintf('============================================================\n');


%% ============================================================
%  RECONSTRUCTION
%  ============================================================

function data = reconstructDehomogenizedField( ...
        meshData,bNodal,rhoNodal,epsMicro,Nx,Ny,nSmoothB,useExactPhases)

    coord = meshData.coord;

    if size(coord,1) ~= numel(bNodal)
        error('b field: %d values for %d mesh nodes.', ...
            numel(bNodal),size(coord,1));
    end

    if size(coord,1) ~= numel(rhoNodal)
        error('rho field: %d values for %d mesh nodes.', ...
            numel(rhoNodal),size(coord,1));
    end

    xMin = min(coord(:,1));
    xMax = max(coord(:,1));
    yMin = min(coord(:,2));
    yMax = max(coord(:,2));

    xGrid = linspace(xMin,xMax,Nx);
    yGrid = linspace(yMin,yMax,Ny);

    [X,Y] = meshgrid(xGrid,yGrid);

    dx = xGrid(2)-xGrid(1);
    dy = yGrid(2)-yGrid(1);

    fprintf('Interpolating optimized fields...\n');

    Fb = scatteredInterpolant( ...
        coord(:,1),coord(:,2),real(bNodal(:)), ...
        'linear','nearest');

    Frho = scatteredInterpolant( ...
        coord(:,1),coord(:,2),real(rhoNodal(:)), ...
        'linear','nearest');

    b   = Fb(X,Y);
    rho = Frho(X,Y);

    b   = real(b);
    rho = real(rho);

    b   = min(max(b,-0.6),0.6);
    rho = min(max(rho,1e-6),1-1e-6);

    if nSmoothB > 0
        b = smoothScalarGrid(b,nSmoothB);
    end

    fprintf('b   min/max = %.6e  %.6e\n', ...
        min(b(:)),max(b(:)));
    fprintf('rho min/max = %.6e  %.6e\n', ...
        min(rho(:)),max(rho(:)));

    if useExactPhases
        fprintf('Using exact phases because b = 0.\n');

        phi1 = (X-xMin)/epsMicro;
        phi2 = (Y-yMin)/epsMicro;

        curl1RMS = 0;
        curl2RMS = 0;
        residual1RMS = 0;
        residual2RMS = 0;
    else
        fprintf('Computing compatible phases with LSQR...\n');

        a = exp(b.^2);
        d = (1+b.^2)./a;

        q1x =  d/epsMicro;
        q1y = -b/epsMicro;

        q2x = -b/epsMicro;
        q2y =  a/epsMicro;

        fprintf('Phase 1:\n');
        phi1 = solveGradientProjectionLSQR(q1x,q1y,dx,dy);

        fprintf('Phase 2:\n');
        phi2 = solveGradientProjectionLSQR(q2x,q2y,dx,dy);

        [dq1x_dy,~] = gradient(q1x,dy,dx);
        [~,dq1y_dx] = gradient(q1y,dy,dx);

        [dq2x_dy,~] = gradient(q2x,dy,dx);
        [~,dq2y_dx] = gradient(q2y,dy,dx);

        curl1 = dq1y_dx-dq1x_dy;
        curl2 = dq2y_dx-dq2x_dy;

        [phi1Y,phi1X] = gradient(phi1,dy,dx);
        [phi2Y,phi2X] = gradient(phi2,dy,dx);

        residual1 = sqrt( ...
            (phi1X-q1x).^2 + ...
            (phi1Y-q1y).^2);

        residual2 = sqrt( ...
            (phi2X-q2x).^2 + ...
            (phi2Y-q2y).^2);

        curl1RMS = sqrt(mean(curl1(:).^2));
        curl2RMS = sqrt(mean(curl2(:).^2));

        residual1RMS = sqrt(mean(residual1(:).^2));
        residual2RMS = sqrt(mean(residual2(:).^2));

        fprintf('RMS curl 1     = %.6e\n',curl1RMS);
        fprintf('RMS curl 2     = %.6e\n',curl2RMS);
        fprintf('RMS residual 1 = %.6e\n',residual1RMS);
        fprintf('RMS residual 2 = %.6e\n',residual2RMS);
    end

    % Coordenadas periodicas centradas em [-0.5,0.5).
    xi1 = phi1-floor(phi1)-0.5;
    xi2 = phi2-floor(phi2)-0.5;

    % Para uma celula quadrada com porosidade definida por rho:
    % area do furo = 1-rho.
    holeHalfSize = sqrt(1-rho)/2;
    holeHalfSize = max(holeHalfSize,1e-6);

    levelSet = max( ...
        abs(xi1)./holeHalfSize, ...
        abs(xi2)./holeHalfSize)-1;

    % levelSet < 0: furo
    % levelSet >= 0: material solido
    solid = double(levelSet >= 0);

    reconstructedVolumeFraction = mean(solid(:));

    fprintf('Grid solid fraction = %.6f\n', ...
        reconstructedVolumeFraction);

    data = struct();

    data.X = X;
    data.Y = Y;

    data.xGrid = xGrid;
    data.yGrid = yGrid;

    data.b   = b;
    data.rho = rho;

    data.phi1 = phi1;
    data.phi2 = phi2;

    data.levelSet = levelSet;
    data.solid    = solid;

    data.epsMicro = epsMicro;
    data.Nx       = Nx;
    data.Ny       = Ny;

    data.reconstructedVolumeFraction = ...
        reconstructedVolumeFraction;

    data.curl1RMS     = curl1RMS;
    data.curl2RMS     = curl2RMS;
    data.residual1RMS = residual1RMS;
    data.residual2RMS = residual2RMS;

end


%% ============================================================
%  LSQR COMPATIBLE-PHASE SOLVER
%  ============================================================

function phi = solveGradientProjectionLSQR(qx,qy,dx,dy)

    [Ny,Nx] = size(qx);

    nNodes = Nx*Ny;

    nHorizontalEdges = Ny*(Nx-1);
    nVerticalEdges   = (Ny-1)*Nx;

    nRows = nHorizontalEdges+nVerticalEdges+1;
    maxNonZeros = 2*(nHorizontalEdges+nVerticalEdges)+1;

    rowIndex = zeros(maxNonZeros,1);
    colIndex = zeros(maxNonZeros,1);
    values   = zeros(maxNonZeros,1);
    rhs      = zeros(nRows,1);

    row = 0;
    nz  = 0;

    for iy = 1:Ny
        for ix = 1:Nx-1

            row = row+1;

            nodeLeft  = sub2ind([Ny,Nx],iy,ix);
            nodeRight = sub2ind([Ny,Nx],iy,ix+1);

            nz = nz+1;
            rowIndex(nz) = row;
            colIndex(nz) = nodeLeft;
            values(nz)   = -1/dx;

            nz = nz+1;
            rowIndex(nz) = row;
            colIndex(nz) = nodeRight;
            values(nz)   = 1/dx;

            rhs(row) = 0.5*(qx(iy,ix)+qx(iy,ix+1));

        end
    end

    for iy = 1:Ny-1
        for ix = 1:Nx

            row = row+1;

            nodeBottom = sub2ind([Ny,Nx],iy,ix);
            nodeTop    = sub2ind([Ny,Nx],iy+1,ix);

            nz = nz+1;
            rowIndex(nz) = row;
            colIndex(nz) = nodeBottom;
            values(nz)   = -1/dy;

            nz = nz+1;
            rowIndex(nz) = row;
            colIndex(nz) = nodeTop;
            values(nz)   = 1/dy;

            rhs(row) = 0.5*(qy(iy,ix)+qy(iy+1,ix));

        end
    end

    anchorWeight = 1e3;

    row = row+1;
    nz  = nz+1;

    rowIndex(nz) = row;
    colIndex(nz) = 1;
    values(nz)   = anchorWeight;

    rhs(row) = 0;

    rowIndex = rowIndex(1:nz);
    colIndex = colIndex(1:nz);
    values   = values(1:nz);

    G = sparse( ...
        rowIndex,colIndex,values,nRows,nNodes);

    fprintf('System: %d x %d, nnz = %d\n', ...
        size(G,1),size(G,2),nnz(G));

    tolLSQR = 1e-6;
    maxIter = 300;

    [phiVector,flag,relres,iter] = ...
        lsqr(G,rhs,tolLSQR,maxIter);

    fprintf('LSQR flag = %d, relres = %.6e, iterations = %d\n', ...
        flag,relres,iter);

    if flag ~= 0 && flag ~= 1
        warning('LSQR returned flag %d.',flag);
    end

    phi = reshape(phiVector,Ny,Nx);

end


%% ============================================================
%  LIGHT SMOOTHING
%  ============================================================

function uSmooth = smoothScalarGrid(u,nSteps)

    uSmooth = u;

    kernel = [
        0 1 0
        1 4 1
        0 1 0
    ];

    kernel = kernel/sum(kernel(:));

    for k = 1:nSteps

        oldValues = uSmooth;

        uSmooth = conv2( ...
            oldValues,kernel,'same');

        uSmooth(1,:)   = oldValues(1,:);
        uSmooth(end,:) = oldValues(end,:);
        uSmooth(:,1)   = oldValues(:,1);
        uSmooth(:,end) = oldValues(:,end);

    end

end


%% ============================================================
%  SAVE MAT, FIGURES AND VTI
%  ============================================================

function saveCaseResults(data,outputFolder,baseName,epsMicro)

    matFile = fullfile(outputFolder,[baseName '.mat']);

    save(matFile,'data','-v7.3');

    % Solid geometry
    figSolid = figure('Color','w');

    imagesc(data.xGrid,data.yGrid,data.solid);
    axis xy;
    axis equal tight;
    box on;

    xlabel('$x$','Interpreter','latex');
    ylabel('$y$','Interpreter','latex');

    title(sprintf( ...
        '%s: dehomogenized solid, $\\varepsilon=%.3f$', ...
        strrep(baseName,'_',' '),epsMicro), ...
        'Interpreter','latex');

    colormap(gray);
    colorbar;

    exportgraphics( ...
        figSolid, ...
        fullfile(outputFolder,[baseName '_Solid.png']), ...
        'Resolution',300);

    savefig( ...
        figSolid, ...
        fullfile(outputFolder,[baseName '_Solid.fig']));

    % ------------------------------------------------------------
    % Preview estilo "material vermelho / vazio azul"
    % ------------------------------------------------------------
    figPreview = figure('Color',[0.35 0.38 0.47]);

    ax = axes(figPreview);
    imagesc(ax,data.xGrid,data.yGrid,data.solid);
    axis(ax,'xy');
    axis(ax,'equal');
    axis(ax,'tight');
    hold(ax,'on');

    % Contorno da fronteira do level set
    contour(ax,data.X,data.Y,data.levelSet,[0 0], ...
        'Color',[0.95 0.55 0.35], ...
        'LineWidth',0.5);

    % Azul para vazio (0), vermelho para sólido (1)
    cmap = [
        0.16 0.28 0.70   % vazio
        0.72 0.00 0.08   % sólido
        ];
    colormap(ax,cmap);
    caxis(ax,[0 1]);

    set(ax,'Color',[0.35 0.38 0.47]);
    set(ax,'XTick',[]);
    set(ax,'YTick',[]);
    box(ax,'off');

    title(ax,'','Interpreter','none');

    exportgraphics( ...
        figPreview, ...
        fullfile(outputFolder,[baseName '_Preview.png']), ...
        'Resolution',300);

    savefig( ...
        figPreview, ...
        fullfile(outputFolder,[baseName '_Preview.fig']));

    % Level-set function
    figLS = figure('Color','w');

    contourf( ...
        data.X,data.Y,data.levelSet,40, ...
        'LineColor','none');

    hold on;

    contour( ...
        data.X,data.Y,data.levelSet,[0 0], ...
        'k','LineWidth',0.8);

    axis equal tight;
    box on;
    colorbar;

    xlabel('$x$','Interpreter','latex');
    ylabel('$y$','Interpreter','latex');

    title(sprintf( ...
        '%s: level set, $\\varepsilon=%.3f$', ...
        strrep(baseName,'_',' '),epsMicro), ...
        'Interpreter','latex');

    exportgraphics( ...
        figLS, ...
        fullfile(outputFolder,[baseName '_LevelSet.png']), ...
        'Resolution',300);

    savefig( ...
        figLS, ...
        fullfile(outputFolder,[baseName '_LevelSet.fig']));

    % VTI for ParaView
    vtiFile = fullfile(outputFolder,[baseName '.vti']);

    fields = struct();

    fields.LevelSet = data.levelSet;
    fields.Solid    = data.solid;
    fields.b        = data.b;
    fields.rho      = data.rho;
    fields.phi1     = data.phi1;
    fields.phi2     = data.phi2;

    writeVTIImageData( ...
        vtiFile,data.xGrid,data.yGrid,fields);

    fprintf('Saved:\n');
    fprintf('  %s\n',matFile);
    fprintf('  %s\n',vtiFile);

end


%% ============================================================
%  VTI WRITER FOR PARAVIEW
%  ============================================================

function writeVTIImageData(fileName,xGrid,yGrid,fields)

    Nx = numel(xGrid);
    Ny = numel(yGrid);

    dx = xGrid(2)-xGrid(1);
    dy = yGrid(2)-yGrid(1);

    x0 = xGrid(1);
    y0 = yGrid(1);

    fid = fopen(fileName,'w');

    if fid < 0
        error('Could not create VTI file: %s',fileName);
    end

    cleanupObject = onCleanup(@() fclose(fid));

    fprintf(fid,'<?xml version="1.0"?>\n');
    fprintf(fid, ...
        '<VTKFile type="ImageData" version="0.1" byte_order="LittleEndian">\n');

    fprintf(fid, ...
        ['  <ImageData WholeExtent="0 %d 0 %d 0 0" ', ...
         'Origin="%.16g %.16g 0" Spacing="%.16g %.16g 1">\n'], ...
        Nx-1,Ny-1,x0,y0,dx,dy);

    fprintf(fid, ...
        '    <Piece Extent="0 %d 0 %d 0 0">\n', ...
        Nx-1,Ny-1);

    fprintf(fid,'      <PointData>\n');

    fieldNames = fieldnames(fields);

    for iField = 1:numel(fieldNames)

        name = fieldNames{iField};
        values = fields.(name);

        if ~isequal(size(values),[Ny,Nx])
            error('VTI field %s has an incompatible size.',name);
        end

        % Transpose so x varies fastest, as expected by VTK.
        valuesVTK = values.';
        valuesVTK = valuesVTK(:);

        fprintf(fid, ...
            '        <DataArray type="Float32" Name="%s" format="ascii">\n', ...
            name);

        chunkSize = 10000;

        for iStart = 1:chunkSize:numel(valuesVTK)

            iEnd = min(iStart+chunkSize-1,numel(valuesVTK));

            fprintf(fid,'%.8g ',valuesVTK(iStart:iEnd));
            fprintf(fid,'\n');

        end

        fprintf(fid,'        </DataArray>\n');

    end

    fprintf(fid,'      </PointData>\n');
    fprintf(fid,'      <CellData/>\n');
    fprintf(fid,'    </Piece>\n');
    fprintf(fid,'  </ImageData>\n');
    fprintf(fid,'</VTKFile>\n');

end
