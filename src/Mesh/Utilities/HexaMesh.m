function mesh = HexaMesh(length, height, width, nx, ny, nz)
    x = linspace(0, length, nx+1);
    y = linspace(0, height, ny+1);
    z = linspace(0, width, nz+1);
    [X, Y, Z] = meshgrid(x, y, z);
    Xr = X(:); Yr = Y(:); Zr = Z(:);
    coor = [Xr, Yr, Zr];
    
    [ix, iy, iz] = ndgrid(1:nx, 1:ny, 1:nz);
    next_x = ny + 1; next_z = (nx+1)*(ny+1);
    
    offset_z = pagetranspose((iz - 1) * next_z);
    offset_x = pagetranspose((ix - 1) * next_x);
    offset_y = pagetranspose(iy - 1);
    
    node1 = offset_z(:) + offset_x(:) + offset_y(:)+1;
    node1 = node1(:);
    node2 = node1 + 1;
    node3 = node1 + next_x + 1;
    node4 = node1 + next_x;
    node5 = node1 + next_z;
    node6 = node2 + next_z;
    node7 = node3 + next_z;
    node8 = node4 + next_z;
    
    glob = [node2,node3,node7,node6,node1,node4,node8,node5];
    
    s.coord = coor;
    s.connec = glob;
    
    mesh = Mesh.create(s);
end