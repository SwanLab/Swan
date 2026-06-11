function printParaview3D(mesh, uFun, thetaFun, wFun, zLayer, outputPath, suffix, varargin)
% printParaview3D - Export 3D extruded shell mesh to VTU for ParaView.
%
% Builds a volumetric mesh by replicating the 2D mid-plane mesh through
% the thickness using FSDT kinematics:
%
%   U_x(x,y,z) = u(x,y) + z * theta_x(x,y)
%   U_y(x,y,z) = v(x,y) + z * theta_y(x,y)
%   U_z(x,y,z) = w(x,y)
%
% The resulting file contains:
%   - 3D wedge elements (triangular prisms) extruded from each triangle
%   - Displacement vector field (Ux, Uy, Uz) at every 3D node
%   - Scalar fields: w, u_x, u_y, theta_x, theta_y
%
% INPUTS:
%   mesh       - Swan mesh object (.coord nNodes x 2, .connec nElem x 3)
%   uFun       - LagrangianFunction with in-plane displacements (nNodes x 2)
%   thetaFun   - LagrangianFunction with rotations (nNodes x 2)
%   wFun       - LagrangianFunction with transverse displacement (nNodes x 1)
%   zLayer     - Cell array of layer interface z-coordinates {z0, z1, ..., zN}
%   outputPath - Folder where the .vtu file will be saved
%   suffix     - String suffix for the filename (e.g., '_Static')
%
% OPTIONAL NAME-VALUE:
%   'nZLayers'      - Number of through-thickness divisions (default: 10)
%   'scaleFactor'   - Amplification factor for displacements (default: 1)
%   'strainFun'     - Cell array of LagrangianFunctions with strain fields
%   'stressFun'     - Cell array of LagrangianFunctions with stress fields
%   'vonMisesFun'   - Cell array of LagrangianFunctions with von Mises
%
% EXAMPLE:
%   printParaview3D(obj.mesh, obj.uFun, obj.thetaFun, obj.wFun, ...
%       obj.zLayer, 'STATIC', '_Static', 'nZLayers', 12);

    % ======================== PARSE INPUTS ========================
    p = inputParser;
    addParameter(p, 'nZLayers', 10, @isnumeric);
    addParameter(p, 'scaleFactor', 1, @isnumeric);
    addParameter(p, 'strainFun', {}, @iscell);
    addParameter(p, 'stressFun', {}, @iscell);
    addParameter(p, 'vonMisesFun', {}, @iscell);
    parse(p, varargin{:});

    nZLayers    = p.Results.nZLayers;
    scaleFactor = p.Results.scaleFactor;
    strainFun   = p.Results.strainFun;
    stressFun   = p.Results.stressFun;
    vonMisesFun = p.Results.vonMisesFun;

    % ======================== 2D MESH DATA ========================
    coord2D  = mesh.coord;          % nNodes2D x 2
    connec2D = mesh.connec;         % nElem x 3 (triangles)
    nNodes2D = size(coord2D, 1);
    nElem    = size(connec2D, 1);

    % Through-thickness coordinates
    zInterfaces = [zLayer{:}];
    zBot = min(zInterfaces);
    zTop = max(zInterfaces);
    zCoords = linspace(zBot, zTop, nZLayers + 1)';
    nZNodes = length(zCoords);  % nZLayers + 1

    % ======================== EXTRACT 2D FIELDS ========================
    u_x     = uFun.fValues(:, 1);       % nNodes2D x 1
    u_y     = uFun.fValues(:, 2);       % nNodes2D x 1
    theta_x = thetaFun.fValues(:, 1);   % nNodes2D x 1
    theta_y = thetaFun.fValues(:, 2);   % nNodes2D x 1
    w_val   = wFun.fValues(:, 1);       % nNodes2D x 1

    % ======================== BUILD 3D NODES ========================
    % Total nodes = nNodes2D * nZNodes
    nNodes3D = nNodes2D * nZNodes;
    coord3D  = zeros(nNodes3D, 3);
    disp3D   = zeros(nNodes3D, 3);

    % Scalar fields for ParaView
    w_3D       = zeros(nNodes3D, 1);
    ux_3D      = zeros(nNodes3D, 1);
    uy_3D      = zeros(nNodes3D, 1);
    thetax_3D  = zeros(nNodes3D, 1);
    thetay_3D  = zeros(nNodes3D, 1);
    zCoord_3D  = zeros(nNodes3D, 1);

    fprintf('  Building 3D mesh: %d nodes x %d z-layers = %d total nodes\n', ...
        nNodes2D, nZNodes, nNodes3D);

    for iZ = 1:nZNodes
        z_k = zCoords(iZ);
        nodeOffset = (iZ - 1) * nNodes2D;
        idx = nodeOffset + (1:nNodes2D);

        % FSDT kinematics: u_3D = {u + z*theta_x, v + z*theta_y, w}
        Ux = u_x + z_k * theta_x;
        Uy = u_y + z_k * theta_y;
        Uz = w_val;

        % Undeformed coordinates (mid-plane x,y extruded to z_k)
        coord3D(idx, 1) = coord2D(:, 1);
        coord3D(idx, 2) = coord2D(:, 2);
        coord3D(idx, 3) = z_k;

        % Displacement vector
        disp3D(idx, 1) = scaleFactor * Ux;
        disp3D(idx, 2) = scaleFactor * Uy;
        disp3D(idx, 3) = scaleFactor * Uz;

        % Scalar fields
        w_3D(idx)      = w_val;
        ux_3D(idx)     = u_x;
        uy_3D(idx)     = u_y;
        thetax_3D(idx) = theta_x;
        thetay_3D(idx) = theta_y;
        zCoord_3D(idx) = z_k;
    end

    % ======================== BUILD 3D CONNECTIVITY ========================
    % Wedge elements (VTK_WEDGE = 13): 6 nodes per element
    % Bottom triangle (layer i) + Top triangle (layer i+1)
    nElem3D  = nElem * nZLayers;
    connec3D = zeros(nElem3D, 6);

    for iZ = 1:nZLayers
        botOffset = (iZ - 1) * nNodes2D;
        topOffset = iZ * nNodes2D;
        elemOffset = (iZ - 1) * nElem;

        for iE = 1:nElem
            n1 = connec2D(iE, 1);
            n2 = connec2D(iE, 2);
            n3 = connec2D(iE, 3);

            % VTK wedge: bottom tri (n1,n2,n3), top tri (n4,n5,n6)
            connec3D(elemOffset + iE, :) = [botOffset + n1, botOffset + n2, botOffset + n3, ...
                                             topOffset + n1, topOffset + n2, topOffset + n3];
        end
    end

    % ======================== STRESS/STRAIN AT 3D NODES ========================
    % If strain/stress functions are provided, interpolate them at each z-level.
    % The user typically provides these evaluated at specific z positions;
    % here we evaluate the FSDT strain field directly for each z-layer.

    hasStrain = ~isempty(strainFun);
    hasStress = ~isempty(stressFun);
    hasVM     = ~isempty(vonMisesFun);

    % ======================== WRITE VTU FILE ========================
    filename = fullfile(outputPath, ['mesh3D' suffix '.vtu']);
    fprintf('  Writing 3D VTU file: %s\n', filename);

    fid = fopen(filename, 'w');
    if fid == -1
        error('Cannot open file: %s', filename);
    end

    % Header
    fprintf(fid, '<?xml version="1.0"?>\n');
    fprintf(fid, '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">\n');
    fprintf(fid, '<UnstructuredGrid>\n');
    fprintf(fid, '<Piece NumberOfPoints="%d" NumberOfCells="%d">\n', nNodes3D, nElem3D);

    % --- Points ---
    fprintf(fid, '<Points>\n');
    fprintf(fid, '<DataArray type="Float64" NumberOfComponents="3" format="ascii">\n');
    fprintf(fid, '%.6e %.6e %.6e\n', coord3D');
    fprintf(fid, '</DataArray>\n');
    fprintf(fid, '</Points>\n');

    % --- Cells ---
    fprintf(fid, '<Cells>\n');

    % Connectivity (0-based)
    fprintf(fid, '<DataArray type="Int32" Name="connectivity" format="ascii">\n');
    fprintf(fid, '%d %d %d %d %d %d\n', (connec3D - 1)');
    fprintf(fid, '</DataArray>\n');

    % Offsets
    fprintf(fid, '<DataArray type="Int32" Name="offsets" format="ascii">\n');
    offsets = (1:nElem3D) * 6;
    fprintf(fid, '%d\n', offsets);
    fprintf(fid, '</DataArray>\n');

    % Types (VTK_WEDGE = 13)
    fprintf(fid, '<DataArray type="UInt8" Name="types" format="ascii">\n');
    fprintf(fid, '%d\n', 13 * ones(nElem3D, 1));
    fprintf(fid, '</DataArray>\n');

    fprintf(fid, '</Cells>\n');

    % --- PointData ---
    fprintf(fid, '<PointData Vectors="Displacement">\n');

    % Displacement vector
    writeDataArray(fid, 'Displacement', disp3D, 3);

    % Scalar fields
    writeDataArray(fid, 'w',       w_3D,      1);
    writeDataArray(fid, 'u_x',     ux_3D,     1);
    writeDataArray(fid, 'u_y',     uy_3D,     1);
    writeDataArray(fid, 'theta_x', thetax_3D, 1);
    writeDataArray(fid, 'theta_y', thetay_3D, 1);
    writeDataArray(fid, 'z',       zCoord_3D, 1);

    % Optional: Strain fields
    if hasStrain && numel(strainFun) >= 2
        nCompStrain = size(strainFun{1}.fValues, 2);
        compNames = {'xx', 'yy', 'yz', 'xz', 'xy'};
        for comp = 1:min(nCompStrain, 5)
            strain_3D = zeros(nNodes3D, 1);
            for iZ = 1:nZNodes
                z_k = zCoords(iZ);
                k = find(zInterfaces <= z_k + 1e-9, 1, 'last');
                if isempty(k), k = 1; end
                if k >= length(zInterfaces), k = length(zInterfaces)-1; end
                
                z0 = zInterfaces(k);
                z1 = zInterfaces(k+1);
                t_interp = (z_k - z0) / (z1 - z0);
                
                idx = (iZ-1)*nNodes2D + (1:nNodes2D);
                strain_3D(idx) = (1-t_interp) * strainFun{k}.fValues(:,comp) ...
                               + t_interp * strainFun{k+1}.fValues(:,comp);
            end
            writeDataArray(fid, sprintf('strain_%s', compNames{comp}), strain_3D, 1);
        end
    end

    % Optional: Stress fields
    if hasStress && numel(stressFun) >= 2
        nCompStress = size(stressFun{1}.fValues, 2);
        compNames = {'xx', 'yy', 'yz', 'xz', 'xy'};
        for comp = 1:min(nCompStress, 5)
            stress_3D = zeros(nNodes3D, 1);
            for iZ = 1:nZNodes
                z_k = zCoords(iZ);
                k = find(zInterfaces <= z_k + 1e-9, 1, 'last');
                if isempty(k), k = 1; end
                if k >= length(zInterfaces), k = length(zInterfaces)-1; end
                
                z0 = zInterfaces(k);
                z1 = zInterfaces(k+1);
                t_interp = (z_k - z0) / (z1 - z0);
                
                idx = (iZ-1)*nNodes2D + (1:nNodes2D);
                stress_3D(idx) = (1-t_interp) * stressFun{k}.fValues(:,comp) ...
                               + t_interp * stressFun{k+1}.fValues(:,comp);
            end
            writeDataArray(fid, sprintf('stress_%s', compNames{comp}), stress_3D, 1);
        end
    end

    % Optional: von Mises at top and bottom surfaces
    if hasVM && numel(vonMisesFun) >= 2
        vm_3D = zeros(nNodes3D, 1);
        for iZ = 1:nZNodes
            z_k = zCoords(iZ);
            k = find(zInterfaces <= z_k + 1e-9, 1, 'last');
            if isempty(k), k = 1; end
            if k >= length(zInterfaces), k = length(zInterfaces)-1; end
            
            z0 = zInterfaces(k);
            z1 = zInterfaces(k+1);
            t_interp = (z_k - z0) / (z1 - z0);
            
            idx = (iZ-1)*nNodes2D + (1:nNodes2D);
            vm_3D(idx) = (1-t_interp) * vonMisesFun{k}.fValues(:,1) ...
                       + t_interp * vonMisesFun{k+1}.fValues(:,1);
        end
        writeDataArray(fid, 'vonMises', vm_3D, 1);
    end

    fprintf(fid, '</PointData>\n');

    % --- CellData (empty) ---
    fprintf(fid, '<CellData>\n');
    fprintf(fid, '</CellData>\n');

    % Footer
    fprintf(fid, '</Piece>\n');
    fprintf(fid, '</UnstructuredGrid>\n');
    fprintf(fid, '</VTKFile>\n');

    fclose(fid);

    fprintf('  ✓ 3D VTU file saved (%d nodes, %d wedge elements)\n', nNodes3D, nElem3D);
    fprintf('  → Open in ParaView and apply "Warp By Vector" (Displacement) to see deformed shape.\n');
    fprintf('  → Use scale factor in Warp to amplify if displacements are small.\n');

end

% ======================== HELPER ========================
function writeDataArray(fid, name, data, nComp)
    fprintf(fid, '<DataArray type="Float64" Name="%s" NumberOfComponents="%d" format="ascii">\n', ...
        name, nComp);
    if nComp == 1
        fprintf(fid, '%.6e\n', data);
    elseif nComp == 3
        fprintf(fid, '%.6e %.6e %.6e\n', data');
    else
        fmt = [repmat('%.6e ', 1, nComp), '\n'];
        fprintf(fid, fmt, data');
    end
    fprintf(fid, '</DataArray>\n');
end
