function m = TriangleMesh(length, height, nx, ny)
    % Generate coordinates
    x1 = linspace(0,length,nx);
    x2 = linspace(0,height,ny);
    % Create the grid
    [xv,yv] = meshgrid(x1,x2);
    % Triangulate the mesh to obtain coordinates and connectivities
    [F,V] = mesh2tri(xv,yv,zeros(size(xv)),'x');
    s.coord  = V(:,1:2);
    s.connec = F;
    m = Mesh.create(s);
end