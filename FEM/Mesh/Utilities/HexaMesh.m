function mesh = HexaMesh(length, height, width, nx, ny, nz)
    % Coord
    x = linspace(0, length, nx+1);
    y = linspace(0, height, ny+1);
    z = linspace(0, width,  nz+1);
    [X,Y,Z] = meshgrid(x,y,z);
    npnod = size(X,1)*size(X,2)*size(X,3);
    Xr = reshape(X, npnod,1);
    Yr = reshape(Y, npnod,1);
    Zr = reshape(Z, npnod,1);
    coor = [Xr, Yr, Zr];
    
    % Connec
    glob = [];
    next_x = ny + 1;
    next_z = (nx+1)*(ny+1);
    for iZ = 1: nz
        z_toadd = (nx+1)*(ny+1) * (iZ-1);
        for iX = 1: nx
            x_toadd = z_toadd + next_x*(iX-1);
            for iY = 1: ny
                y_line = x_toadd + [iY, iY + 1];
                z0_plane = [y_line, flip(next_x + y_line) ];
                nods = [z0_plane, next_z + z0_plane];
                connec = [nods(2), nods(3), nods(7), nods(6), ...
                    nods(1), nods(4), nods(8), nods(5)];
                glob = [glob; connec];
            end
        end
    end
    
    s.coord = coor;
    s.connec = glob;
    
    mesh = Mesh.create(s);
end