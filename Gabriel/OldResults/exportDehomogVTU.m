function exportDehomogVTU(mesh, lsValues, fileName)
% EXPORTDEHOMOGVTU Exporta a level set da dehomogeneizacao para VTU.
%
% Uso:
%   exportDehomogVTU(mesh,lsValues,'Dehomogenization_eps005.vtu')
%
% O arquivo contém:
%   - level_set: valor médio da level set por elemento;
%   - solid_positive: 1 quando level_set >= 0;
%   - solid_negative: 1 quando level_set <= 0.
%
% No ParaView:
%   1. Abra o arquivo .vtu;
%   2. Apply;
%   3. Color By -> level_set;
%   4. Para extrair uma fase, use Threshold;
%   5. Para visualizar a interface, use Contour com valor 0.

    arguments
        mesh
        lsValues (:,1) double
        fileName (1,:) char
    end

    %% ------------------------------------------------------------
    % Coordenadas
    %% ------------------------------------------------------------

    coord = mesh.coord;

    if size(coord,2) == 2
        points = [coord, zeros(size(coord,1),1)];
    elseif size(coord,2) == 3
        points = coord;
    else
        error('mesh.coord deve possuir duas ou três colunas.');
    end

    nPoints = size(points,1);

    %% ------------------------------------------------------------
    % Conectividade
    %% ------------------------------------------------------------

    connec = mesh.connec;

    % Em algumas classes, connec pode estar transposta.
    if size(connec,1) <= 4 && size(connec,2) > size(connec,1)
        connec = connec.';
    end

    nCells      = size(connec,1);
    nNodesCell  = size(connec,2);

    if min(connec(:)) == 0
        connecZero = connec;
    else
        connecZero = connec - 1;
    end

    if max(connecZero(:)) >= nPoints
        error('A conectividade contém índices maiores que o número de nós.');
    end

    %% ------------------------------------------------------------
    % Tipo de célula VTK
    %% ------------------------------------------------------------

    switch nNodesCell
        case 3
            vtkCellType = 5;  % VTK_TRIANGLE
        case 4
            vtkCellType = 9;  % VTK_QUAD
        otherwise
            error(['Tipo de elemento não suportado: %d nós por elemento. ', ...
                   'Esta função aceita triângulos ou quadriláteros.'], ...
                   nNodesCell);
    end

    offsets = (nNodesCell:nNodesCell:nNodesCell*nCells).';
    types   = vtkCellType*ones(nCells,1);

    %% ------------------------------------------------------------
    % Converter lsValues para um valor por elemento
    %% ------------------------------------------------------------

    lsValues = real(lsValues(:));

    if numel(lsValues) == nCells*nNodesCell

        % Caso P1D: um valor por nó local de cada elemento.
        lsByElement = reshape(lsValues,nNodesCell,nCells);
        levelSetCell = mean(lsByElement,1).';

    elseif numel(lsValues) == nCells

        % Caso P0: um valor por elemento.
        levelSetCell = lsValues;

    elseif numel(lsValues) == nPoints

        % Caso P1: um valor global por nó.
        levelSetCell = mean( ...
            reshape(lsValues(connec(:)),nCells,nNodesCell),2);

    else
        error([ ...
            'Não consegui associar lsValues à malha.\n', ...
            'numel(lsValues) = %d\n', ...
            'nPoints         = %d\n', ...
            'nCells          = %d\n', ...
            'nCells*nNodesCell = %d'], ...
            numel(lsValues),nPoints,nCells,nCells*nNodesCell);
    end

    % As duas convenções são exportadas para evitar dúvida sobre o sinal.
    solidPositive = double(levelSetCell >= 0);
    solidNegative = double(levelSetCell <= 0);

    %% ------------------------------------------------------------
    % Escrever VTU ASCII
    %% ------------------------------------------------------------

    fid = fopen(fileName,'w');

    if fid < 0
        error('Não foi possível criar o arquivo: %s',fileName);
    end

    cleaner = onCleanup(@() fclose(fid));

    fprintf(fid,'<?xml version="1.0"?>\n');
    fprintf(fid,'<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">\n');
    fprintf(fid,'  <UnstructuredGrid>\n');
    fprintf(fid,'    <Piece NumberOfPoints="%d" NumberOfCells="%d">\n', ...
        nPoints,nCells);

    %% Points

    fprintf(fid,'      <Points>\n');
    fprintf(fid,'        <DataArray type="Float64" NumberOfComponents="3" format="ascii">\n');

    fprintf(fid,'          %.16e %.16e %.16e\n',points.');

    fprintf(fid,'        </DataArray>\n');
    fprintf(fid,'      </Points>\n');

    %% Cells

    fprintf(fid,'      <Cells>\n');

    fprintf(fid,'        <DataArray type="Int32" Name="connectivity" format="ascii">\n');
    for iCell = 1:nCells
        fprintf(fid,'          ');
        fprintf(fid,'%d ',connecZero(iCell,:));
        fprintf(fid,'\n');
    end
    fprintf(fid,'        </DataArray>\n');

    fprintf(fid,'        <DataArray type="Int32" Name="offsets" format="ascii">\n');
    fprintf(fid,'          %d\n',offsets);
    fprintf(fid,'        </DataArray>\n');

    fprintf(fid,'        <DataArray type="UInt8" Name="types" format="ascii">\n');
    fprintf(fid,'          %d\n',types);
    fprintf(fid,'        </DataArray>\n');

    fprintf(fid,'      </Cells>\n');

    %% Cell data

    fprintf(fid,'      <CellData Scalars="level_set">\n');

    fprintf(fid,'        <DataArray type="Float64" Name="level_set" format="ascii">\n');
    fprintf(fid,'          %.16e\n',levelSetCell);
    fprintf(fid,'        </DataArray>\n');

    fprintf(fid,'        <DataArray type="UInt8" Name="solid_positive" format="ascii">\n');
    fprintf(fid,'          %d\n',solidPositive);
    fprintf(fid,'        </DataArray>\n');

    fprintf(fid,'        <DataArray type="UInt8" Name="solid_negative" format="ascii">\n');
    fprintf(fid,'          %d\n',solidNegative);
    fprintf(fid,'        </DataArray>\n');

    fprintf(fid,'      </CellData>\n');

    fprintf(fid,'    </Piece>\n');
    fprintf(fid,'  </UnstructuredGrid>\n');
    fprintf(fid,'</VTKFile>\n');

    fprintf('\nVTU salvo com sucesso:\n%s\n',fullfile(pwd,fileName));
    fprintf('Pontos:    %d\n',nPoints);
    fprintf('Elementos: %d\n',nCells);
    fprintf('Level set min/max: %.6e %.6e\n', ...
        min(levelSetCell),max(levelSetCell));
end