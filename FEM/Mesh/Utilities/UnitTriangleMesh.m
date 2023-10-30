function m = UnitTriangleMesh(nx, ny)
    % Generate coordinates
    x1 = linspace(0,1,nx);
    x2 = linspace(0,1,ny);
    % Create the grid
    [xv,yv] = meshgrid(x1,x2);
    % Triangulate the mesh to obtain coordinates and connectivities
    [F,V] = mesh2tri(xv,yv,zeros(size(xv)),'x');
    s.coord  = V(:,1:2);
    s.connec = F;
    m = Mesh(s);
end