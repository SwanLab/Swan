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
                vert_line = x_toadd + [iY, iY + 1];
                plane_elem = [vert_line, flip(next_x + vert_line) ];
                connec = [plane_elem, next_z + plane_elem];
                glob = [glob; connec];
            end
        end
    end
    
    s.coord = coor;
    s.connec = glob;
    
    mesh = Mesh.create(s);
end