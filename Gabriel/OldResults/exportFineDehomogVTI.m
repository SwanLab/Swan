function exportFineDehomogVTI(phaseFile, outputFile, NxFine, NyFine)
% EXPORTFINEDEHOMOGVTI
%
% Reconstrói a level set da microestrutura em uma grade regular fina
% independente da malha da topology e exporta para ParaView.
%
% Exemplo:
%
% exportFineDehomogVTI( ...
%   'CompatiblePhaseResults_eps005_300x200.mat', ...
%   'DehomogenizationFine_eps005.vti', ...
%   1200,600);
%
% Campos exportados:
%   level_set
%   solid_positive
%   solid_negative
%   rho
%   phi1_periodic
%   phi2_periodic

    arguments
        phaseFile  (1,:) char
        outputFile (1,:) char
        NxFine     (1,1) double {mustBeInteger,mustBeGreaterThan(NxFine,1)}
        NyFine     (1,1) double {mustBeInteger,mustBeGreaterThan(NyFine,1)}
    end

    %% ------------------------------------------------------------
    % Carregar fases
    %% ------------------------------------------------------------

    if ~exist(phaseFile,'file')
        error('Arquivo não encontrado: %s',phaseFile);
    end

    P = load(phaseFile);

    if ~isfield(P,'phaseResults')
        error('O arquivo não contém phaseResults.');
    end

    R = P.phaseResults;

    requiredFields = {'X','Y','phi1','phi2','rho'};

    for i = 1:numel(requiredFields)
        if ~isfield(R,requiredFields{i})
            error('Campo ausente em phaseResults: %s',requiredFields{i});
        end
    end

    %% ------------------------------------------------------------
    % Grade original das fases
    %% ------------------------------------------------------------

    xPhase = R.X(1,:);
    yPhase = R.Y(:,1);

    xMin = min(xPhase);
    xMax = max(xPhase);
    yMin = min(yPhase);
    yMax = max(yPhase);

    %% ------------------------------------------------------------
    % Interpoladores
    %% ------------------------------------------------------------

    Fphi1 = griddedInterpolant( ...
        {yPhase,xPhase},R.phi1, ...
        'linear','nearest');

    Fphi2 = griddedInterpolant( ...
        {yPhase,xPhase},R.phi2, ...
        'linear','nearest');

    Frho = griddedInterpolant( ...
        {yPhase,xPhase},R.rho, ...
        'linear','nearest');

    %% ------------------------------------------------------------
    % Grade fina independente
    %% ------------------------------------------------------------

    xFine = linspace(xMin,xMax,NxFine);
    yFine = linspace(yMin,yMax,NyFine);

    [Xfine,Yfine] = meshgrid(xFine,yFine);

    fprintf('\n============================================================\n');
    fprintf('FINE DEHOMOGENIZATION EXPORT\n');
    fprintf('============================================================\n');
    fprintf('Phase file   : %s\n',phaseFile);
    fprintf('Fine grid    : %d x %d\n',NxFine,NyFine);
    fprintf('Fine points  : %d\n',NxFine*NyFine);

    fprintf('Interpolating phase fields...\n');

    phi1 = Fphi1(Yfine,Xfine);
    phi2 = Fphi2(Yfine,Xfine);
    rho  = Frho(Yfine,Xfine);

    rho = min(max(real(rho),1e-6),1-1e-6);

    %% ------------------------------------------------------------
    % Coordenadas periódicas locais
    %% ------------------------------------------------------------

    xi1 = phi1-floor(phi1)-0.5;
    xi2 = phi2-floor(phi2)-0.5;

    %% ------------------------------------------------------------
    % Level set do furo quadrado
    %
    % rho = fração sólida
    % lado total do furo = sqrt(1-rho)
    % metade do lado = sqrt(1-rho)/2
    %% ------------------------------------------------------------

    holeHalfSize = sqrt(1-rho)/2;
    holeHalfSize = max(holeHalfSize,1e-8);

    levelSet = max( ...
        abs(xi1)./holeHalfSize, ...
        abs(xi2)./holeHalfSize)-1;

    solidPositive = double(levelSet >= 0);
    solidNegative = double(levelSet <= 0);

    fprintf('Level set min/max: %.6e %.6e\n', ...
        min(levelSet(:)),max(levelSet(:)));

    %% ------------------------------------------------------------
    % Escrever VTI ASCII
    %% ------------------------------------------------------------

    dx = (xMax-xMin)/(NxFine-1);
    dy = (yMax-yMin)/(NyFine-1);

    fid = fopen(outputFile,'w');

    if fid < 0
        error('Não foi possível criar: %s',outputFile);
    end

    cleaner = onCleanup(@() fclose(fid));

    fprintf(fid,'<?xml version="1.0"?>\n');
    fprintf(fid,['<VTKFile type="ImageData" version="0.1" ', ...
                 'byte_order="LittleEndian">\n']);

    fprintf(fid,[ ...
        '  <ImageData WholeExtent="0 %d 0 %d 0 0" ', ...
        'Origin="%.16e %.16e 0" ', ...
        'Spacing="%.16e %.16e 1">\n'], ...
        NxFine-1,NyFine-1,xMin,yMin,dx,dy);

    fprintf(fid,'    <Piece Extent="0 %d 0 %d 0 0">\n', ...
        NxFine-1,NyFine-1);

    fprintf(fid,'      <PointData Scalars="level_set">\n');

    writeDataArray(fid,'level_set',levelSet);
    writeDataArray(fid,'rho',rho);
    writeDataArray(fid,'phi1_periodic',xi1);
    writeDataArray(fid,'phi2_periodic',xi2);
    writeDataArray(fid,'solid_positive',solidPositive);
    writeDataArray(fid,'solid_negative',solidNegative);

    fprintf(fid,'      </PointData>\n');
    fprintf(fid,'      <CellData></CellData>\n');
    fprintf(fid,'    </Piece>\n');
    fprintf(fid,'  </ImageData>\n');
    fprintf(fid,'</VTKFile>\n');

    fprintf('\nArquivo VTI salvo:\n%s\n',fullfile(pwd,outputFile));
    fprintf('============================================================\n');

end


function writeDataArray(fid,name,A)

    fprintf(fid,[ ...
        '        <DataArray type="Float64" ', ...
        'Name="%s" NumberOfComponents="1" format="ascii">\n'], ...
        name);

    % VTK espera x variando mais rapidamente.
    % A do MATLAB é Ny x Nx; A(:) já percorre primeiro y.
    % A transposta coloca x como índice mais rápido.
    values = A.';

    fprintf(fid,'          %.16e\n',values(:));

    fprintf(fid,'        </DataArray>\n');

end